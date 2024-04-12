import re

def convert_files(input_file_path,file_name,cutoff=400,rel_cutoff=60):
  input_file = input_file_path
  output_file_path = input_file[:-4] + '.inp'

  project_name = file_name



  output_file = open(output_file_path, 'w')
  output_file.write(f"""&GLOBAL
  PROJECT_NAME {project_name}
  RUN_TYPE GEO_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL
&MOTION
&GEO_OPT
  OPTIMIZER BFGS
  MAX_ITER 1000
  MAX_DR  [bohr] 0.003
  &BFGS
  &END
&END GEO_OPT
&END MOTION
&FORCE_EVAL
  &DFT
    BASIS_SET_FILE_NAME /rds/general/user/ewolpert/home/cp2k/basis_set/BASIS_MOLOPT
    POTENTIAL_FILE_NAME /rds/general/user/ewolpert/home/cp2k/basis_set/GTH_POTENTIALS
  &QS
   EPS_DEFAULT 1.0E-10
  &END QS
  &MGRID
    NGRIDS 4
    CUTOFF {cutoff}
    REL_CUTOFF {rel_cutoff}
  &END MGRID
  &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 300
      &OT
        PRECONDITIONER FULL_SINGLE_INVERSE
        MINIMIZER DIIS
      &END OT
      &OUTER_SCF ! repeat the inner SCF cycle 10 times
        MAX_SCF 10
        EPS_SCF 1.0E-6 ! must match the above
      &END      
      &PRINT
        &RESTART OFF
        &END
      &END
  &END SCF
 &XC
  &VDW_POTENTIAL
    POTENTIAL_TYPE PAIR_POTENTIAL
    &PAIR_POTENTIAL
      TYPE DFTD3
      CALCULATE_C9_TERM .TRUE.
      PARAMETER_FILE_NAME /rds/general/user/ewolpert/home/cp2k/basis_set/dftd3.dat
      REFERENCE_FUNCTIONAL PBE
      R_CUTOFF 10.
      EPS_CN 0.01
    &END PAIR_POTENTIAL
  &END VDW_POTENTIAL
  &XC_FUNCTIONAL PBE
  &END XC_FUNCTIONAL
&END XC
&END DFT
&SUBSYS
 &CELL
  ABC 60 60 60 
  ALPHA_BETA_GAMMA 90 90 90
 &END CELL
&COORD  
""")
  input_file = f"{input_file_path}"
  with open(input_file, 'r') as input:
      lines = input.readlines()
  for line in lines[2:]:
      output_file.write(line)
  output_file.write("""&END COORD
&KIND C
  BASIS_SET TZVP-MOLOPT-GTH
  POTENTIAL GTH-PBE-q4
&END KIND
&KIND O
  BASIS_SET TZVP-MOLOPT-GTH
  POTENTIAL GTH-PBE-q6
&END KIND
&KIND F
  BASIS_SET TZVP-MOLOPT-GTH
  POTENTIAL GTH-PBE-q7
&END KIND
&KIND S
  BASIS_SET TZVP-MOLOPT-GTH
  POTENTIAL GTH-PBE-q6
&END KIND
&KIND Br
  BASIS_SET DZVP-MOLOPT-SR-GTH
  POTENTIAL GTH-PBE-q7
&END KIND
&KIND Cl
  BASIS_SET DZVP-MOLOPT-SR-GTH
  POTENTIAL GTH-PBE-q7
&END KIND
&KIND N
  BASIS_SET TZVP-MOLOPT-GTH
  POTENTIAL GTH-PBE-q5
&END KIND
&KIND H
  BASIS_SET TZVP-MOLOPT-GTH
  POTENTIAL GTH-PBE-q1
&END KIND
&END SUBSYS
&END FORCE_EVAL""")
pass

if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    convert_files(input_file, output_file)


