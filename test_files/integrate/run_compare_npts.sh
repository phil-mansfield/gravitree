go run simulate.go --method="brute" --gen="circ" --npts="3" --dt="1e-4" --saveto="vary_npts";
go run simulate.go --method="brute" --gen="circ" --npts="4" --dt="1e-4" --saveto="vary_npts";
go run simulate.go --method="brute" --gen="circ" --npts="5" --dt="1e-4" --saveto="vary_npts";

go run simulate.go --method="brute" --gen="circ" --npts="3" --dt="1e-4" --saveto="vary_npts" --int="rk4";
go run simulate.go --method="brute" --gen="circ" --npts="4" --dt="1e-4" --saveto="vary_npts" --int="rk4";
go run simulate.go --method="brute" --gen="circ" --npts="5" --dt="1e-4" --saveto="vary_npts" --int="rk4";