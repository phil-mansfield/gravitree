# Varies the eta in adaptive time stepping.

rm -r vary_eta

etas=("0.025" "0.0025" "0.00025" "0.000025")

for eta in "${etas[@]}"; do
    go run simulate.go --method="brute" --gen="dipole" --int="lfadp" --feps="1" --saveto="vary_eta" --dt="1e-2" --dp=0. --eta="$eta"
    go run simulate.go --method="brute" --gen="dipole" --int="lfadp" --feps="1" --saveto="vary_eta" --dt="1e-2" --dp=.1 --eta="$eta"
    go run simulate.go --method="brute" --gen="dipole" --int="lfadp" --feps="1" --saveto="vary_eta" --dt="1e-2" --dp=.5 --eta="$eta"
    go run simulate.go --method="brute" --gen="dipole" --int="lfadp" --feps="1" --saveto="vary_eta" --dt="1e-2" --dp=.9 --eta="$eta"
done
