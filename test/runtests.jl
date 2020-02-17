using SciBase, Test, Unitful

println("Running tests:")


x0 = 0
x = 2
@test x_abs(x0,x) == 2

println("All tests passed!")