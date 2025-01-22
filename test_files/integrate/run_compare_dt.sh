# time resolution vs energy error
go run simulate.go --method="brute" --gen="single" --dt="1e-1" --saveto="vary_dt"
go run simulate.go --method="brute" --gen="single" --int="rk4" --dt="1e-1" --saveto="vary_dt"

go run simulate.go --method="brute" --gen="single" --dt="1e-2" --saveto="vary_dt"
go run simulate.go --method="brute" --gen="single" --int="rk4" --dt="1e-2" --saveto="vary_dt"

go run simulate.go --method="brute" --gen="single" --dt="1e-3" --saveto="vary_dt"
go run simulate.go --method="brute" --gen="single" --int="rk4" --dt="1e-3" --saveto="vary_dt"

go run simulate.go --method="brute" --gen="single" --dt="1e-4" --saveto="vary_dt"
go run simulate.go --method="brute" --gen="single" --int="rk4" --dt="1e-4" --saveto="vary_dt"