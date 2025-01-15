# go run simulate_tracers.go --method="approx" --gen="stream"
# go run simulate_tracers.go --method="brute" --gen="stream"

go run simulate_tracers.go --method="approx" --gen="circ" --npts="3"
go run simulate_tracers.go --method="brute" --gen="circ" --npts="3"

python3 plot_circ.py 3

go run simulate_tracers.go --method="approx" --gen="circ" --npts="4"
go run simulate_tracers.go --method="brute" --gen="circ" --npts="4"

python3 plot_circ.py 4

go run simulate_tracers.go --method="approx" --gen="circ" --npts="5"
go run simulate_tracers.go --method="brute" --gen="circ" --npts="5"

python3 plot_circ.py 5

# go run simulate_tracers.go --method="approx" --gen="sing"
# go run simulate_tracers.go --method="brute" --gen="sing"
# python3 plot_sing.py