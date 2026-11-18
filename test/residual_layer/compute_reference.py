#!/usr/bin/env python3
"""Reference solutions for the residual layer and gas dry deposition test case.

Writes ref_rl_gas.txt, ref_rl_aero.txt (with the residual layer) and
ref_norl_gas.txt, ref_norl_aero.txt (without it), in the column layout of
extract_gas and extract_aero_time. The gases are passive (no chemistry) and the
aerosol has no coagulation, so every quantity follows from the input files by
mixing, emission and first-order deposition alone. The time stepping is the
operator splitting of PartMC: in each step of 60 s ending at time t, emissions
and deposition use the profile values at t and the mixing height h(t), and
the mixed layer exchanges air with the layers above when h changes from
h(t - 60 s) to h(t). Without a residual layer, a growing mixed layer entrains
background air; with it, the air between the old and new mixing height goes
to the residual layer when the mixed layer becomes shallower and comes back
before any background air when it grows again.
"""
import math

DT, T_MAX, T_OUT = 60.0, 32400.0, 1800.0
P, TEMP, R = 1e5, 288.0, 8.314472                     # PartMC's gas constant
NAIR = P / (R * TEMP)                                # mol m-3
HEIGHT = [(0, 1000), (3600, 1000), (7200, 200), (14400, 200), (21600, 1500), (32400, 1500)]
GAS = ['H2SO4', 'HNO3', 'HCl', 'NH3', 'NO', 'NO2']
GAS_INIT = dict(NO=10.0, HNO3=4.0, HCl=3.0)
GAS_BACK = dict(NO2=5.0)
GAS_EMIT = dict(NH3=1e-8)                            # mol m-2 s-1 at rate 1
EMIT_RATE = [(0, 0), (7200, 1), (14400, 0)]          # gas and aerosol emissions
DEPO = [(0, 1, 0.01), (3600, 1, 0.0), (7200, 1, 0.005), (14400, 1, 0.0), (25200, 2, 0.01)]   # HNO3: (time, rate, v_d)
# aerosol modes: mono-disperse, one species each (aero_data.dat: density)
MODES = dict(SI=(3e-5, 10000.0), SB=(4e-5, 100.0), SE=(5e-5, 1500.0))
AERO = ['SI', 'SB', 'SE']
NUM_INIT = dict(SI=1e9)                              # m-3
NUM_BACK = dict(SB=2e9)                              # m-3
NUM_EMIT = dict(SE=1e8)                              # m-2 s-1 at rate 1

def height(t):
    for (t0, h0), (t1, h1) in zip(HEIGHT[:-1], HEIGHT[1:]):
        if t0 <= t <= t1:
            return h0 + (h1 - h0) * (t - t0) / (t1 - t0)
    return HEIGHT[-1][1]

def step_value(profile, t):
    """Value of a step profile at t: the last set-point at or before t."""
    v = profile[0][1:]
    for p in profile:
        if p[0] <= t:
            v = p[1:]
    return v

def particle_mass(sp):
    d, rho = MODES[sp]
    return math.pi / 6 * d ** 3 * rho

def run(residual_layer):
    c = {g: GAS_INIT.get(g, 0.0) for g in GAS}                         # ppb
    n = {a: NUM_INIT.get(a, 0.0) for a in AERO}                        # m-3
    bg_c = {g: GAS_BACK.get(g, 0.0) for g in GAS}
    bg_n = {a: NUM_BACK.get(a, 0.0) for a in AERO}
    rl_c = {g: 0.0 for g in GAS}; rl_n = {a: 0.0 for a in AERO}; top = height(0)
    rows = [(0.0, dict(c), dict(n))]
    n_steps = int(round(T_MAX / DT))
    for i in range(1, n_steps + 1):
        t = i * DT; h0 = height(t - DT); h1 = height(t)
        rate = step_value(EMIT_RATE, t)[0]
        for g, e in GAS_EMIT.items():
            c[g] += rate * e * DT / h1 * 1e9 / NAIR
        for a, e in NUM_EMIT.items():
            n[a] += rate * e * DT / h1
        if not residual_layer and h1 > h0:
            f = h0 / h1
            c = {g: f * c[g] + (1 - f) * bg_c[g] for g in GAS}
            n = {a: f * n[a] + (1 - f) * bg_n[a] for a in AERO}
        r, v = step_value(DEPO, t)
        c['HNO3'] *= math.exp(-r * v * DT / h1)
        if residual_layer:
            if h1 < h0:
                d = max(top - h0, 0.0); f = (h0 - h1) / (d + h0 - h1)
                rl_c = {g: (1 - f) * rl_c[g] + f * c[g] for g in GAS}
                rl_n = {a: (1 - f) * rl_n[a] + f * n[a] for a in AERO}
                top = max(top, h0)
            elif h1 > h0:
                d_rl = max(min(h1, top) - h0, 0.0); d_bg = h1 - h0 - d_rl
                c = {g: (h0 * c[g] + d_rl * rl_c[g] + d_bg * bg_c[g]) / h1 for g in GAS}
                n = {a: (h0 * n[a] + d_rl * rl_n[a] + d_bg * bg_n[a]) / h1 for a in AERO}
                if h1 >= top:
                    top = h1; rl_c = {g: 0.0 for g in GAS}; rl_n = {a: 0.0 for a in AERO}
        if abs(t / T_OUT - round(t / T_OUT)) < 1e-9:
            rows.append((t, dict(c), dict(n)))
    return rows

def write(rows, prefix):
    with open(prefix + '_gas.txt', 'w') as f:
        for t, c, n in rows:
            f.write(''.join('%30.15e' % x for x in [t] + [c[g] for g in GAS]) + '\n')
    with open(prefix + '_aero.txt', 'w') as f:
        for t, c, n in rows:
            m = [n[a] * particle_mass(a) for a in AERO]
            f.write(''.join('%30.15e' % x for x in [t, sum(n.values()), sum(m)] + m) + '\n')

if __name__ == '__main__':
    write(run(True), 'ref_rl')
    write(run(False), 'ref_norl')
