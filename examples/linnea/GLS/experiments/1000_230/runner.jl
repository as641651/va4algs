using MatrixGenerator
using LinearAlgebra.BLAS
BLAS.set_num_threads(4)

include("experiments/algorithm36.jl")
include("experiments/algorithm67.jl")
include("experiments/algorithm94.jl")
include("experiments/algorithm53.jl")
include("experiments/algorithm12.jl")
include("experiments/algorithm43.jl")
include("experiments/algorithm84.jl")
include("experiments/algorithm77.jl")
include("experiments/algorithm26.jl")
include("experiments/algorithm57.jl")
include("experiments/algorithm63.jl")
include("experiments/algorithm32.jl")
include("experiments/algorithm90.jl")
include("experiments/algorithm80.jl")
include("experiments/algorithm22.jl")
include("experiments/algorithm73.jl")
include("experiments/algorithm47.jl")
include("experiments/algorithm16.jl")
include("experiments/algorithm56.jl")
include("experiments/algorithm62.jl")
include("experiments/algorithm33.jl")
include("experiments/algorithm91.jl")
include("experiments/algorithm81.jl")
include("experiments/algorithm23.jl")
include("experiments/algorithm72.jl")
include("experiments/algorithm46.jl")
include("experiments/algorithm17.jl")
include("experiments/algorithm37.jl")
include("experiments/algorithm66.jl")
include("experiments/algorithm95.jl")
include("experiments/algorithm52.jl")
include("experiments/algorithm13.jl")
include("experiments/algorithm42.jl")
include("experiments/algorithm85.jl")
include("experiments/algorithm76.jl")
include("experiments/algorithm27.jl")
include("experiments/algorithm59.jl")
include("experiments/algorithm18.jl")
include("experiments/algorithm49.jl")
include("experiments/algorithm1.jl")
include("experiments/algorithm69.jl")
include("experiments/algorithm38.jl")
include("experiments/algorithm28.jl")
include("experiments/algorithm79.jl")
include("experiments/algorithm5.jl")
include("experiments/algorithm68.jl")
include("experiments/algorithm39.jl")
include("experiments/algorithm29.jl")
include("experiments/algorithm78.jl")
include("experiments/algorithm4.jl")
include("experiments/algorithm58.jl")
include("experiments/algorithm19.jl")
include("experiments/algorithm48.jl")
include("experiments/algorithm0.jl")
include("experiments/algorithm98.jl")
include("experiments/algorithm88.jl")
include("experiments/algorithm7.jl")
include("experiments/algorithm3.jl")
include("experiments/algorithm2.jl")
include("experiments/algorithm99.jl")
include("experiments/algorithm89.jl")
include("experiments/algorithm6.jl")
include("experiments/algorithm55.jl")
include("experiments/algorithm92.jl")
include("experiments/algorithm30.jl")
include("experiments/algorithm61.jl")
include("experiments/algorithm71.jl")
include("experiments/algorithm20.jl")
include("experiments/algorithm82.jl")
include("experiments/algorithm14.jl")
include("experiments/algorithm45.jl")
include("experiments/algorithm96.jl")
include("experiments/algorithm65.jl")
include("experiments/algorithm34.jl")
include("experiments/algorithm51.jl")
include("experiments/algorithm41.jl")
include("experiments/algorithm10.jl")
include("experiments/algorithm24.jl")
include("experiments/algorithm75.jl")
include("experiments/algorithm86.jl")
include("experiments/algorithm9.jl")
include("experiments/algorithm97.jl")
include("experiments/algorithm64.jl")
include("experiments/algorithm35.jl")
include("experiments/algorithm50.jl")
include("experiments/algorithm40.jl")
include("experiments/algorithm11.jl")
include("experiments/algorithm25.jl")
include("experiments/algorithm74.jl")
include("experiments/algorithm87.jl")
include("experiments/algorithm8.jl")
include("experiments/algorithm54.jl")
include("experiments/algorithm93.jl")
include("experiments/algorithm31.jl")
include("experiments/algorithm60.jl")
include("experiments/algorithm70.jl")
include("experiments/algorithm21.jl")
include("experiments/algorithm83.jl")
include("experiments/algorithm15.jl")
include("experiments/algorithm44.jl")

include("operand_generator.jl")

function main()

    matrices = operand_generator()

    io = open("/Users/aravind/exercise/performance-analyzer/VariantsCompare/dev-notes/linnea/GLS/experiments/1000_230/run_times.csv","w")
    write(io, "case:concept:name;concept:name;concept:flops;concept:operation;concept:kernel;timestamp:start;timestamp:end\n")

    n = 2000
    rand(n, n)*rand(n, n) # this seems to help to reduce some startup noise

    ret,times = algorithm36(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm36_to_eventlog(io, "algorithm36", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm67(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm67_to_eventlog(io, "algorithm67", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm94(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm94_to_eventlog(io, "algorithm94", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm53(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm53_to_eventlog(io, "algorithm53", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm12(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm12_to_eventlog(io, "algorithm12", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm43(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm43_to_eventlog(io, "algorithm43", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm84(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm84_to_eventlog(io, "algorithm84", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm77(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm77_to_eventlog(io, "algorithm77", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm26(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm26_to_eventlog(io, "algorithm26", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm57(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm57_to_eventlog(io, "algorithm57", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm63(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm63_to_eventlog(io, "algorithm63", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm32(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm32_to_eventlog(io, "algorithm32", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm90(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm90_to_eventlog(io, "algorithm90", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm80(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm80_to_eventlog(io, "algorithm80", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm22(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm22_to_eventlog(io, "algorithm22", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm73(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm73_to_eventlog(io, "algorithm73", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm47(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm47_to_eventlog(io, "algorithm47", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm16(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm16_to_eventlog(io, "algorithm16", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm56(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm56_to_eventlog(io, "algorithm56", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm62(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm62_to_eventlog(io, "algorithm62", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm33(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm33_to_eventlog(io, "algorithm33", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm91(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm91_to_eventlog(io, "algorithm91", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm81(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm81_to_eventlog(io, "algorithm81", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm23(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm23_to_eventlog(io, "algorithm23", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm72(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm72_to_eventlog(io, "algorithm72", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm46(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm46_to_eventlog(io, "algorithm46", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm17(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm17_to_eventlog(io, "algorithm17", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm37(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm37_to_eventlog(io, "algorithm37", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm66(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm66_to_eventlog(io, "algorithm66", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm95(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm95_to_eventlog(io, "algorithm95", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm52(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm52_to_eventlog(io, "algorithm52", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm13(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm13_to_eventlog(io, "algorithm13", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm42(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm42_to_eventlog(io, "algorithm42", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm85(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm85_to_eventlog(io, "algorithm85", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm76(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm76_to_eventlog(io, "algorithm76", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm27(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm27_to_eventlog(io, "algorithm27", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm59(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm59_to_eventlog(io, "algorithm59", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm18(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm18_to_eventlog(io, "algorithm18", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm49(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm49_to_eventlog(io, "algorithm49", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm1(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm1_to_eventlog(io, "algorithm1", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm69(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm69_to_eventlog(io, "algorithm69", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm38(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm38_to_eventlog(io, "algorithm38", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm28(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm28_to_eventlog(io, "algorithm28", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm79(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm79_to_eventlog(io, "algorithm79", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm5(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm5_to_eventlog(io, "algorithm5", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm68(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm68_to_eventlog(io, "algorithm68", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm39(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm39_to_eventlog(io, "algorithm39", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm29(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm29_to_eventlog(io, "algorithm29", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm78(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm78_to_eventlog(io, "algorithm78", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm4(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm4_to_eventlog(io, "algorithm4", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm58(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm58_to_eventlog(io, "algorithm58", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm19(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm19_to_eventlog(io, "algorithm19", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm48(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm48_to_eventlog(io, "algorithm48", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm0(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm0_to_eventlog(io, "algorithm0", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm98(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm98_to_eventlog(io, "algorithm98", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm88(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm88_to_eventlog(io, "algorithm88", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm7(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm7_to_eventlog(io, "algorithm7", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm3(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm3_to_eventlog(io, "algorithm3", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm2(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm2_to_eventlog(io, "algorithm2", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm99(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm99_to_eventlog(io, "algorithm99", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm89(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm89_to_eventlog(io, "algorithm89", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm6(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm6_to_eventlog(io, "algorithm6", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm55(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm55_to_eventlog(io, "algorithm55", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm92(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm92_to_eventlog(io, "algorithm92", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm30(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm30_to_eventlog(io, "algorithm30", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm61(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm61_to_eventlog(io, "algorithm61", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm71(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm71_to_eventlog(io, "algorithm71", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm20(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm20_to_eventlog(io, "algorithm20", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm82(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm82_to_eventlog(io, "algorithm82", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm14(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm14_to_eventlog(io, "algorithm14", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm45(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm45_to_eventlog(io, "algorithm45", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm96(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm96_to_eventlog(io, "algorithm96", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm65(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm65_to_eventlog(io, "algorithm65", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm34(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm34_to_eventlog(io, "algorithm34", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm51(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm51_to_eventlog(io, "algorithm51", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm41(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm41_to_eventlog(io, "algorithm41", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm10(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm10_to_eventlog(io, "algorithm10", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm24(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm24_to_eventlog(io, "algorithm24", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm75(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm75_to_eventlog(io, "algorithm75", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm86(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm86_to_eventlog(io, "algorithm86", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm9(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm9_to_eventlog(io, "algorithm9", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm97(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm97_to_eventlog(io, "algorithm97", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm64(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm64_to_eventlog(io, "algorithm64", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm35(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm35_to_eventlog(io, "algorithm35", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm50(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm50_to_eventlog(io, "algorithm50", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm40(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm40_to_eventlog(io, "algorithm40", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm11(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm11_to_eventlog(io, "algorithm11", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm25(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm25_to_eventlog(io, "algorithm25", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm74(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm74_to_eventlog(io, "algorithm74", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm87(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm87_to_eventlog(io, "algorithm87", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm8(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm8_to_eventlog(io, "algorithm8", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm54(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm54_to_eventlog(io, "algorithm54", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm93(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm93_to_eventlog(io, "algorithm93", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm31(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm31_to_eventlog(io, "algorithm31", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm60(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm60_to_eventlog(io, "algorithm60", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm70(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm70_to_eventlog(io, "algorithm70", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm21(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm21_to_eventlog(io, "algorithm21", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm83(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm83_to_eventlog(io, "algorithm83", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm15(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm15_to_eventlog(io, "algorithm15", times)
    temp = rand(25000) # cache trashing

    ret,times = algorithm44(map(MatrixGenerator.unwrap, map(copy, matrices))...)
    write_algorithm44_to_eventlog(io, "algorithm44", times)
    temp = rand(25000) # cache trashing




end

main()