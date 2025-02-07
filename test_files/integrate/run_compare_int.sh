# comparing integrators

feps=("0" "0.01" "0.1" "0.5" "0.75" "1.0")

for feps_value in "${feps[@]}"; do
    go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="$feps_value" --saveto="vary_feps" --dt="1e-3"
done

# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="0" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-2" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-1" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="10" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="100" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1000" --saveto="vary_feps" --dt="1e-3"

# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="0" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1e-2" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1e-1" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="10" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="100" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1000" --saveto="vary_feps" --dt="1e-2"

# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="0" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-2" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-1" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="10" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="100" --saveto="vary_feps" --dt="1e-2"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1000" --saveto="vary_feps" --dt="1e-2"

# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="0" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1e-2" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1e-1" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="10" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="100" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1000" --saveto="vary_feps" --dt="1e-4"

# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="0" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-2" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-1" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="10" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="100" --saveto="vary_feps" --dt="1e-4"
# go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1000" --saveto="vary_feps" --dt="1e-4"

# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="leapfrog" --feps="0" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="leapfrog" --feps="1e-2" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="leapfrog" --feps="1e-1" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="leapfrog" --feps="1" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="leapfrog" --feps="10" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="leapfrog" --feps="100" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="leapfrog" --feps="1000" --saveto="vary_feps" --dt="1e-3"

# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="rk4" --feps="0" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="rk4" --feps="1e-2" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="rk4" --feps="1e-1" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="rk4" --feps="1" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="rk4" --feps="10" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="rk4" --feps="100" --saveto="vary_feps" --dt="1e-3"
# go run simulate.go --method="brute" --gen="ics" --npts="4" --int="rk4" --feps="1000" --saveto="vary_feps" --dt="1e-3"