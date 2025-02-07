go run simulate.go --method="brute" --gen="circ" --npts="3" --int="leapfrog" --dt="1e-2" --saveto="vary_theta" --theta="1e-3" --crit="pkd"
go run simulate.go --method="brute" --gen="circ" --npts="4" --int="leapfrog" --dt="1e-2" --saveto="vary_theta" --theta="1e-2" --crit="pkd"
go run simulate.go --method="brute" --gen="circ" --npts="5" --int="leapfrog" --dt="1e-2" --saveto="vary_theta" --theta="1e-1" --crit="pkd"

go run simulate.go --method="approx" --gen="circ" --npts="3" --int="leapfrog" --dt="1e-2" --saveto="vary_theta" --theta="1e-3" --crit="pkd"
go run simulate.go --method="approx" --gen="circ" --npts="4" --int="leapfrog" --dt="1e-2" --saveto="vary_theta" --theta="1e-2" --crit="pkd"
go run simulate.go --method="approx" --gen="circ" --npts="5" --int="leapfrog" --dt="1e-2" --saveto="vary_theta" --theta="1e-1" --crit="pkd"