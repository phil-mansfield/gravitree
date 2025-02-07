# generates dipole simulations

rm -r vary_dipole

dt_values=("1e-1" "1e-2" "1e-3" "1e-4" "1e-5" "1e-6")
# dp_values=(0.0 0.025 0.05 0.075 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
dp_values=($(seq 0 0.025 0.5))


for dt in "${dt_values[@]}"; do
    for dp in "${dp_values[@]}"; do
        go run simulate.go --method="brute" --gen="dipole" --int="leapfrog" --feps="0" --saveto="vary_dipole" --dt="$dt" --dp=$dp
    done
done