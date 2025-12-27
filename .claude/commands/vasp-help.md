Get help with VASP parameters and settings.

Arguments: $ARGUMENTS (parameter name or category)

Based on $ARGUMENTS, provide detailed information:

1. If a specific parameter is given (e.g., "encut", "ispin", "ibrion"):
   - Explain what it does
   - Show valid values
   - Give usage examples
   - Show how to set it in this interface

2. If a category is given:
   - "relaxation" - IBRION, ISIF, NSW, EDIFFG, POTIM
   - "electronic" - ENCUT, ISMEAR, SIGMA, EDIFF, NELM
   - "magnetic" - ISPIN, MAGMOM, LORBIT
   - "hybrid" - LHFCALC, HFSCREEN, AEXX, ALGO
   - "vdw" - IVDW, VDW_S6, VDW_SR
   - "dft+u" - LDAU, LDAUTYPE, LDAUL, LDAUU, LDAUJ
   - "neb" - IMAGES, SPRING, LCLIMB, ICHAIN
   - "md" - MDALGO, TEBEG, TEEND, SMASS, POTIM

3. If no argument, show the main categories and ask what the user needs help with.

Also show how to use the parameter preset functions:
- get_vdw_params('d3bj')
- get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=4.0)})
- get_hybrid_params('hse06')
- get_soc_params()
