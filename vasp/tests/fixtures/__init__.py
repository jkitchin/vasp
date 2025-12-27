"""Test fixtures for VASP calculator tests."""

import os

FIXTURES_DIR = os.path.dirname(__file__)


def get_fixture_path(filename: str) -> str:
    """Get absolute path to a fixture file."""
    return os.path.join(FIXTURES_DIR, filename)


# Sample OUTCAR content for a completed Si calculation
MOCK_OUTCAR_COMPLETE = """
 vasp.6.3.0 18Jan22 (build Apr 01 2022 16:49:06) complex

 executed on             LinuxIFC date 2024.01.15  10:30:45
 running on   32 total cores
 distrk:  each k-point on   32 cores,    1 groups
 distr:  one band on NCORES_PER_BAND=   8 cores,    4 groups


--------------------------------------------------------------------------------------------------------


 POTCAR:    PAW_PBE Si 05Jan2001

   VRHFIN =Si: s2p2
   LEXCH  = PE
   EATOM  =   103.0669 eV,    7.5752 Ry

   TITEL  = PAW_PBE Si 05Jan2001
   LULTRA =        F    use ultrasoft PP ?

 POSCAR = Si2
   PREC = accurate
   ISTART =      0    job   : new run
   ICHARG =      2    charge: atomic
   ISPIN  =      1    non spin polarized
   LNONCOLLINEAR =      F non collinear calculations
   LSORBIT =      F    spin-orbit coupling

 Electronic Coverage (with calculation)
   NELECT =       8.0000    total number of electrons
   NBAND  =        8    number of bands

   ENCUT  =  400.0 eV

 k-points           NKPTS =     56   k-points in BZ     NKDIM =     56   number of bands

 Orbital magnetization related:
   ORBITALMAG=     F  switch on orbital magnetization

 exchange correlation table for  LEXCH =        8
   RHO(1)=    0.500       N(1)  =     2000
   RHO(2)=  100.500       N(2)  =     4000

 POTLOK:  cpu time      0.0154: real time      0.0158

 energy without entropy =      -10.84532612  energy(sigma->0) =      -10.84016421

  free energy    TOTEN  =       -10.84274516 eV

  energy  without entropy=      -10.84532612  energy(sigma->0) =      -10.84016421

 FORCES acting on ions:
    electron-ion (+-)
 ---------------------------------------
   0.15234  -0.12341   0.00123
  -0.15234   0.12341  -0.00123
 ---------------------------------------

  POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000         0.012345     -0.023456      0.001234
      1.35750      1.35750      1.35750        -0.012345      0.023456     -0.001234
 -----------------------------------------------------------------------------------

  STRESS after sorting (KINETIC / LATTICE / TOTAL) in kB:
  -----------------------------------------------------------------------------------
      1.234     -0.123      0.012
     -0.123      1.234      0.012
      0.012      0.012      1.234
  -----------------------------------------------------------------------------------

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  -----------------------------------------------------------------------------------
  Total:      0.00456     0.00456     0.00456    -0.00012     0.00001     0.00001

  in kB       1.23400     1.23400     1.23400    -0.12300     0.01200     0.01200

 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-loss (STOPCAR, key: LABORT)
  energy-cutoff  :      400.00
  volume of cell :       40.89

 number of electron       8.0000000 magnetization
 augmentation part        0.8974212 magnetization

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =         3.42871879
  Ewald energy   TEWEN  =      -115.34682178
  -Hartree energ DENC   =       -27.37234567
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =        13.23958412
  PAW double counting   =        85.45678932      -84.12345678
  entropy T*S    EENTRO =        -0.00516190
  eigenvalues    EBANDS =       -18.95678912
  atomic energy  EATOM  =       206.13368224
  ---------------------------------------------------
  free energy    TOTEN  =       -10.84274516 eV

  energy without entropy =      -10.84532612  energy(sigma->0) =      -10.84016421

 magnetization (x)

 eigenvalue-Loss (STOPCAR, key: LABORT)
   1   0.0000   0.0000   0.0000          1       0.000000
   2   0.0000   0.0000   0.0000          1       0.000000

 E-fermi :   5.4321     XC(G=0):  -9.8765     alpha+bet : -5.4321

 General timing and accounting informance:
 =============================
   Total CPU time used (sec):       45.123
   User time (sec):                 42.456
   System time (sec):                2.667
            Elapsed time (sec):     47.890

           Job finished at 10:31:33  15Jan2024

 General timing and accounting
"""


MOCK_OUTCAR_SPIN = """
 vasp.6.3.0 18Jan22 (build Apr 01 2022 16:49:06) complex

 ISPIN  =      2    spin polarized calculation

 number of electron       8.0000000 magnetization       2.0000000
 augmentation part        0.8974212 magnetization       0.5234123

 magnetization (x)

# PAW_PBE Fe 06Sep2000
   spin component 1
     ion      1     2.34567
   spin component 2
     ion      1     0.12345

  free energy    TOTEN  =      -16.78901234 eV

 E-fermi :   6.1234     XC(G=0): -10.5678     alpha+bet : -6.7890

 General timing and accounting informance:
   Total CPU time used (sec):       95.678
 General timing and accounting
"""


MOCK_VASPRUN_XML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
  <generator>
    <i name="program" type="string">vasp</i>
    <i name="version" type="string">6.3.0</i>
    <i name="subversion" type="string">18Jan22</i>
  </generator>
  <incar>
    <i type="string" name="PREC">accurate</i>
    <i type="int" name="ENCUT">400</i>
    <i type="int" name="ISMEAR">1</i>
    <v type="int" name="KPOINT">8 8 8</v>
    <i type="logical" name="LWAVE">F</i>
    <i type="logical" name="LCHARG">F</i>
  </incar>
  <atominfo>
    <atoms>2</atoms>
    <types>1</types>
    <array name="atoms">
      <dimension dim="1">ion</dimension>
      <field>element</field>
      <field>atomtype</field>
      <set>
        <rc><c>Si</c><c>   1</c></rc>
        <rc><c>Si</c><c>   1</c></rc>
      </set>
    </array>
  </atominfo>
  <structure name="finalpos">
    <crystal>
      <varray name="basis">
        <v>   5.43000000   0.00000000   0.00000000</v>
        <v>   0.00000000   5.43000000   0.00000000</v>
        <v>   0.00000000   0.00000000   5.43000000</v>
      </varray>
      <i name="volume">160.103007</i>
    </crystal>
    <varray name="positions">
      <v>   0.00000000   0.00000000   0.00000000</v>
      <v>   0.25000000   0.25000000   0.25000000</v>
    </varray>
  </structure>
  <calculation>
    <scstep>
      <energy>
        <i name="e_fr_energy">-10.84274516</i>
        <i name="e_wo_entrp">-10.84532612</i>
        <i name="e_0_energy">-10.84016421</i>
      </energy>
    </scstep>
    <energy>
      <i name="e_fr_energy">-10.84274516</i>
      <i name="e_wo_entrp">-10.84532612</i>
      <i name="e_0_energy">-10.84016421</i>
    </energy>
    <varray name="forces">
      <v>   0.01234500  -0.02345600   0.00123400</v>
      <v>  -0.01234500   0.02345600  -0.00123400</v>
    </varray>
    <varray name="stress">
      <v>   1.23400000  -0.12300000   0.01200000</v>
      <v>  -0.12300000   1.23400000   0.01200000</v>
      <v>   0.01200000   0.01200000   1.23400000</v>
    </varray>
    <dos>
      <i name="efermi">5.4321</i>
      <total>
        <array>
          <dimension dim="1">gridpoints</dimension>
          <dimension dim="2">spin</dimension>
          <field>energy</field>
          <field>total</field>
          <field>integrated</field>
          <set>
            <set comment="spin 1">
              <r>  -5.0000   0.0000   0.0000</r>
              <r>  -4.0000   0.1234   0.1234</r>
              <r>  -3.0000   0.5678   0.6912</r>
              <r>  -2.0000   1.2345   1.9257</r>
              <r>  -1.0000   2.3456   4.2713</r>
              <r>   0.0000   3.4567   7.7280</r>
              <r>   1.0000   2.3456  10.0736</r>
              <r>   2.0000   1.2345  11.3081</r>
              <r>   3.0000   0.5678  11.8759</r>
              <r>   4.0000   0.1234  11.9993</r>
              <r>   5.0000   0.0000  12.0000</r>
            </set>
          </set>
        </array>
      </total>
    </dos>
    <eigenvalues>
      <array>
        <dimension dim="1">kpoint</dimension>
        <dimension dim="2">band</dimension>
        <dimension dim="3">spin</dimension>
        <field>eigene</field>
        <field>occ</field>
        <set>
          <set comment="spin 1">
            <set comment="kpoint 1">
              <r>  -5.6789   1.0000</r>
              <r>   2.3456   1.0000</r>
              <r>   2.3457   1.0000</r>
              <r>   2.3458   1.0000</r>
              <r>   6.7890   0.0000</r>
              <r>   6.7891   0.0000</r>
              <r>   6.7892   0.0000</r>
              <r>   8.9012   0.0000</r>
            </set>
            <set comment="kpoint 2">
              <r>  -4.5678   1.0000</r>
              <r>   1.2345   1.0000</r>
              <r>   3.4567   1.0000</r>
              <r>   3.4568   1.0000</r>
              <r>   5.6789   0.0000</r>
              <r>   7.8901   0.0000</r>
              <r>   7.8902   0.0000</r>
              <r>   9.0123   0.0000</r>
            </set>
          </set>
        </set>
      </array>
    </eigenvalues>
    <projected>
      <kpointlist>
        <varray name="kpointlist">
          <v>   0.00000000   0.00000000   0.00000000</v>
          <v>   0.12500000   0.00000000   0.00000000</v>
        </varray>
        <varray name="weights">
          <v>   0.01562500</v>
          <v>   0.09375000</v>
        </varray>
      </kpointlist>
    </projected>
  </calculation>
</modeling>
"""


MOCK_INCAR = """PREC = accurate
ENCUT = 400
ISMEAR = 1
SIGMA = 0.1
EDIFF = 1E-6
LWAVE = .FALSE.
LCHARG = .FALSE.
NSW = 0
IBRION = -1
"""


MOCK_POSCAR = """Si2 diamond
1.0
   5.4300000000   0.0000000000   0.0000000000
   0.0000000000   5.4300000000   0.0000000000
   0.0000000000   0.0000000000   5.4300000000
Si
2
Direct
   0.0000000000   0.0000000000   0.0000000000
   0.2500000000   0.2500000000   0.2500000000
"""


MOCK_KPOINTS = """Automatic mesh
0
Gamma
  8  8  8
  0  0  0
"""


MOCK_POTCAR_HEADER = """  PAW_PBE Si 05Jan2001
   VRHFIN =Si: s2p2
   LEXCH  = PE
   EATOM  =   103.0669 eV,    7.5752 Ry

   TITEL  = PAW_PBE Si 05Jan2001
   LULTRA =        F    use ultrasoft PP ?
   IUNSCR =        1    unscreen: 0-lin 1-nonlin 2-no
   RPACOR =    1.500    partial core radius
   POMASS =   28.085; ZVAL   =    4.000    mass and valenz
   RCORE  =    1.900    outmost cutoff radius
   RWIGS  =    2.480; RWIGS  =    1.312    wigner-seitz radius (au A)
   ENMAX  =  245.345; ENMIN  =  184.009 eV
   End of Dataset
"""


MOCK_OUTCAR_VIBRATIONS = """
 vasp.6.3.0 18Jan22 (build Apr 01 2022 16:49:06) complex

 IBRION =      5    ionic relaxation

 Eigenvectors and calculation

   1 f  =    15.123456 THz   94.98765 2PiTHz  504.56789 cm-1    62.54321 meV
             X         Y         Z           dx          dy          dz
      0.00000   0.00000   0.00000     0.123456   -0.234567    0.012345
      1.35750   1.35750   1.35750    -0.123456    0.234567   -0.012345

   2 f  =    14.567890 THz   91.23456 2PiTHz  486.12345 cm-1    60.23456 meV
             X         Y         Z           dx          dy          dz
      0.00000   0.00000   0.00000     0.234567    0.123456    0.345678
      1.35750   1.35750   1.35750    -0.234567   -0.123456   -0.345678

   3 f/i=    12.345678 THz   77.56789 2PiTHz  411.78901 cm-1    51.01234 meV
             X         Y         Z           dx          dy          dz
      0.00000   0.00000   0.00000     0.111111    0.222222    0.333333
      1.35750   1.35750   1.35750    -0.111111   -0.222222   -0.333333

  free energy    TOTEN  =       -10.84274516 eV

 General timing and accounting informance:
   Total CPU time used (sec):      125.678
 General timing and accounting
"""


# Expected parsed values for assertions
EXPECTED_ENERGY = -10.84274516
EXPECTED_FORCES = [
    [0.012345, -0.023456, 0.001234],
    [-0.012345, 0.023456, -0.001234]
]
EXPECTED_FERMI = 5.4321
