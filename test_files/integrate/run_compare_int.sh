# comparing integrators

go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="0" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1e-2" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1e-1" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="1" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --feps="10" --saveto="vary_feps" --dt="1e-3"

go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="0" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-2" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1e-1" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="1" --saveto="vary_feps" --dt="1e-3"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="rk4" --feps="10" --saveto="vary_feps" --dt="1e-3"
