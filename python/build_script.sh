# On x86
go build -buildmode=c-shared -o gravitree.so gravitree_wrapper.go
# On ARM
go build -buildmode=c-shared -o gravitree.so gravitree_wrapper.go
