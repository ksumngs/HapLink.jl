@testitem "Linkage disequilibrium" begin
    using HapLink: linkage, _correlation_coeff

    genotypes = [474 142; 611 773]

    # https://pbgworks.org/sites/pbgworks.org/files/measuresoflinkagedisequilibrium-111119214123-phpapp01_0.pdf
    @test trunc(linkage(genotypes); digits=4) ≈ 0.0699
    @test trunc(_correlation_coeff(genotypes); digits=2) ≈ 0.30
    @test trunc(significance(genotypes); digits=4) <= 0.0001
end #testitem
