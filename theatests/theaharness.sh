#
# Method to print out general and script specific options
#
print_usage() {

cat >&2 <<EOF
Usage: $0 [options]

OPTIONS
  -c ................ Cleanup (remove generated files)
  -h ................ help: print this message
  -i <input file> ... Namelist file to use
  -I <vars=vals> .... String of VAR=VAL format to modify namelist; e.g.,
                      -I 'mx=45 nonlinear=.true.'
  -n <integer> ...... Override the number of processors to use
  -v ................ verbose
  -V ................ Run valgrind on tests
EOF

  if declare -f extrausage > /dev/null; then extrausage; fi
  exit $1
}
###
##  Arguments for overriding things
#
diff_flags=""
while getopts "chi:I:n:vV" arg
do
  case $arg in
    c ) cleanup=true         ;;  
    h ) print_usage; exit    ;;  
    i ) input_file="-i $OPTARG" ;;  
    I ) input_varstr="-I $OPTARG" ;;  
    n ) nprocs="$OPTARG"      ;;  
    v ) verbose=true         ;;  
    V ) valgrind_cmd="valgrind -q --tool=memcheck --leak-check=yes --num-callers=20 --track-origins=yes" 
        mpiexeccmd="$mpiexeccmd $nprocs $valgrind" ;;
    *)  # To take care of any extra args
      if test -n "$OPTARG"; then
        eval $arg=\"$OPTARG\"
      else
        eval $arg=found
      fi
      ;;
  esac
done
shift $(( $OPTIND - 1 ))
