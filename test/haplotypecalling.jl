@testitem "Linkage disequilibrium" begin
    using HapLink: linkage

    # https://pbgworks.org/sites/pbgworks.org/files/measuresoflinkagedisequilibrium-111119214123-phpapp01_0.pdf
    @test trunc(linkage([474 142; 611 773]); digits=4) ≈ 0.0699
end #testitem
