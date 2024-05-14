import Aqua
import SamApp2024
using Test: @testset

@testset verbose = true "aqua" begin
    Aqua.test_all(SamApp2024; ambiguities = false)
end
