#!/bin/sh

cat <<ENDINFO

Barrel case
---------------------

This simulates a barrel study with aerosol
coagulation.

ENDINFO
sleep 1

echo ../../build/partmc barrel_with_coag.spec
../../build/partmc barrel_with_coag.spec
