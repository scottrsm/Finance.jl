using Finance
import PropCheck
import Random
using Test

# Set random seed.
Random.seed!(1)

# TEST PARAMETERS
N = 20          # Length of data series.
TOL  = 1.0e-10  # Floating point tolerance.
TOL2 = 1.0e-10  # Floating point tolerance.
RNDS = rand(N)  # Data series.

# Sample Standard deviation of data series.
std_gold1 = 0.3281428080051744e0

# Dictionary of EMA values of the RNDS series using window lengths based on the key.
ema_gold = Dict{Int,Vector{Float64}}()
ema_gold[4] = [0.0491718221481211, 0.07647220156972791, 0.2028552023589405, 0.14804900498517393, 0.37007283482556547, 0.5497896629068475, 0.3844389620446481, 0.60366558615181, 0.6874996801680648, 0.7590692711631398, 0.582205431178262, 0.4639364589644358, 0.3235322488226564, 0.14500651100061113, 0.3202461312911789, 0.5053297826661038, 0.6513597594808751, 0.7180832041793237, 0.5856723479410886, 0.6398109909864496]
ema_gold[5] = [0.0491718221481211, 0.07404391268445873, 0.1891855152731569, 0.13925416882623337, 0.34152962850882335, 0.5114791931083412, 0.3852245480494335, 0.5521144166573246, 0.6878872765955792, 0.7598207656882588, 0.538180864704961, 0.49878425282784744, 0.3695047122894986, 0.2113399627735, 0.3040560642390054, 0.49131235212233054, 0.6110810685176296, 0.654940979403851, 0.5862387991746807, 0.65415318677649]
ema_gold[6] = [0.0491718221481211, 0.06840085700911294, 0.15908389615346957, 0.12951112315473623, 0.2897180655700219, 0.4376861018227501, 0.3680106596801032, 0.5192804326092977, 0.6162703580180132, 0.7324976466240967, 0.571828719230447, 0.4967318033109158, 0.44286255318, 0.294910490641581, 0.3391306050324787, 0.42825604120252536, 0.5511587940906428, 0.5954779847157526, 0.5545687083059226, 0.6536399647444989]
ema_gold[7] = [0.0491718221481211, 0.06716362392786011, 0.1520119419461167, 0.12434193803189089, 0.2742408576870402, 0.4126883459380966, 0.3474959557319084, 0.49353069101721414, 0.6019221519797247, 0.6869176136880574, 0.5795515881354328, 0.5141547114199814, 0.4199819538224159, 0.3309939039267013, 0.3713822212175776, 0.45801505782391005, 0.5245899699632272, 0.5795378435432132, 0.5316599082434561, 0.6121109784366969]
ema_gold[8] = [0.0491718221481211, 0.06400175999546284, 0.13463881892199814, 0.11572081736333462, 0.2414707155396819, 0.36326388637904017, 0.32136919260592356, 0.4491607324080409, 0.5533673503954608, 0.6479755303428287, 0.553346577636453, 0.5377975618264139, 0.4576535681719037, 0.3458733652885238, 0.41181065743040823, 0.4708418487535503, 0.5216135610775359, 0.5316949501689585, 0.5134904065388175, 0.5784251823097781]
ema_gold[9] = [0.0491718221481211, 0.0632548720857186, 0.1303343998995003, 0.11236917553964487, 0.23178586587882521, 0.34744510360762515, 0.3076603741906429, 0.42901588189145856, 0.5279742875201132, 0.6213384287441491, 0.5452846088361767, 0.5119256565262541, 0.46944889248205957, 0.3671088770967002, 0.3954647929829185, 0.49022555158084535, 0.5376679701823374, 0.5497792181744611, 0.4945908207325262, 0.566806905159303]
ema_gold[10] = [0.0491718221481211, 0.061237716852839234, 0.1190670572462412, 0.10569079876302581, 0.20930147770647886, 0.3125588041053841, 0.28503629879411285, 0.3937129851628374, 0.48567497699777856, 0.5744290467623708, 0.5218056559899439, 0.5091390979263087, 0.4605263175876842, 0.4012867244928605, 0.4216726617624846, 0.4722306296745527, 0.5458570637517266, 0.5524336358627324, 0.5046395446114016, 0.5351855792411268]

# Dictionary of Moving Std values of the RNDS series using window lengths based on the key.
ema_std_gold = Dict{Int,Vector{Float64}}()
ema_std_gold[4] = [0.2004876151136763, 0.17150857475998477, 0.20156870029598967, 0.1926991864789593, 0.28790232440243957, 0.29032618342019423, 0.320854306484565, 0.3247151074624288, 0.2705236669278982, 0.23399310770564666, 0.3666700827509274, 0.30621327004538423, 0.26768311963119906, 0.24193405493771677, 0.23417697342160243, 0.2911368371206138, 0.26514438184409017, 0.2211498440374814, 0.23559106876707492, 0.1943083385369941]
ema_std_gold[5] = [0.3276599916072654, 0.27725767125845446, 0.27201287204201013, 0.24201489378048516, 0.3150735236451221, 0.3179603342187367, 0.3364122720853384, 0.34641929288650464, 0.3077576771188212, 0.2546163331534271, 0.3402729605104603, 0.2919403558524953, 0.25761717348636465, 0.25385754528339805, 0.2885057991356663, 0.30125243528431817, 0.28205424148780706, 0.23202828811764856, 0.2479817825334222, 0.20933235348688045]
ema_std_gold[6] = [0.36640602984475856, 0.32776660497376964, 0.32238350385541087, 0.2937206309165007, 0.35143343300335034, 0.36748737779047724, 0.36607474665013756, 0.3803473791040501, 0.3567160132026443, 0.3295064740055784, 0.3693199074447034, 0.3264319125989771, 0.3133460000986224, 0.3102810430291659, 0.30644178812968165, 0.34597224646538, 0.3170431054054486, 0.28172899967782206, 0.2692712339169174, 0.2362854549900182]
ema_std_gold[7] = [0.34685040672822903, 0.31036108852350236, 0.3074794587361891, 0.27957622942971805, 0.340947060614906, 0.36248635911678917, 0.35403339265911393, 0.37450567238594995, 0.3589261685016636, 0.33403168416950213, 0.3853250943452138, 0.3401425469383805, 0.31513119683395435, 0.3271537684066131, 0.30330307751985847, 0.32570411645731145, 0.32504040488015506, 0.2718965052647351, 0.2598232112259467, 0.2431565907075897]
ema_std_gold[8] = [0.3796194610752972, 0.3491907517390243, 0.3451261439269724, 0.31975241028546386, 0.36911795520549817, 0.39379253008649767, 0.37945230456176837, 0.402542499879443, 0.3958229466882905, 0.3823218194080913, 0.4018430037199607, 0.379761239558348, 0.3540434930164778, 0.3514794853877449, 0.3294495242640888, 0.3287120783875891, 0.3245277609352351, 0.29846341700232776, 0.26656670309311775, 0.25961552367917917]
ema_std_gold[9] = [0.39029213899400744, 0.35890898611224265, 0.35266758034205065, 0.32618318271941643, 0.3723141907457161, 0.3969116117467207, 0.3791494639323386, 0.4039325066818461, 0.39983157884697323, 0.3889560495899872, 0.40686628638866756, 0.3769443948062868, 0.36897775307110053, 0.36614637750176854, 0.3362917780840278, 0.33935536471103833, 0.3228697162013628, 0.2928978674753235, 0.27132853451453415, 0.2533558393785408]
ema_std_gold[10] = [0.39657167182108954, 0.37089544160218185, 0.3666404749231677, 0.343981039104234, 0.38527608853630707, 0.41183287491297227, 0.39405955033758333, 0.42050672943078116, 0.42247295103177757, 0.4182140373512705, 0.4254030966151893, 0.4030644884019153, 0.388788280153648, 0.4017604591524959, 0.3671809594003703, 0.3582859695943349, 0.35028582672602776, 0.31420261800605814, 0.29090126465225447, 0.2833816102767438]

# Pi to the 11th power.
Random.seed!(1)
RNDS2 = rand(1000)
pi_11_gold = 294204.0179738905e0
r_1000_l1_gold = 0.987213065673524e0
r_1000_l9_gold = 0.06972438782723225e0

#------------------------------------------------------------------------
#--------------      TESTS       ----------------------------------------
#------------------------------------------------------------------------
#
@testset "Finance (Fidelity)            " begin

    @test length(detect_ambiguities(Finance)) == 0
end

@testset "Finance (Moving Average Stats)" begin
    window_sizes = [w for w in 4:10]
    for w in window_sizes
        ema_res = ema(RNDS, w)
        ema_std_res = ema_std(RNDS, w)
        @test length(ema_res) == length(RNDS)
        @test length(ema_std_res) == length(RNDS)
        @test ema_res ≈ ema_gold[w] rtol = TOL
        @test ema_std_res ≈ ema_std_gold[w] rtol = TOL
        @test std(RNDS) ≈ std_gold1 rtol = TOL

        # Compute ema stats.
        stats = ema_stats(RNDS, w)
        @test size(stats) == (length(RNDS), 4)
        @test stats[:, 1] ≈ ema_res rtol = TOL
        @test stats[:, 2] ≈ ema_std_res rtol = TOL

    end
end

@testset "Finance (Entropy Index)       " begin
    # Compute fast exponential
    @test pow_n(2, 7) == 128
    @test pow_n(2, 5) == 32
    @test pow_n(2, 8) == 256
    @test pow_n(π, 11) ≈ pi_11_gold rtol = TOL

    # Check entropy_index
    @test entropy_index(RNDS2; n=100, λ=1.0) ≈ r_1000_l1_gold rtol = TOL
    @test entropy_index(RNDS2; n=100, λ=0.9) ≈ r_1000_l9_gold rtol = TOL
end


@testset "Finance (Test Input Contracts)" begin
    @test_throws DomainError ema(RNDS, 1)
    @test_throws DomainError ema(RNDS, 2)
    @test_throws DomainError ema(RNDS, 3)

    @test_throws DomainError ema_std(RNDS, 1)
    @test_throws DomainError ema_std(RNDS, 2)
    @test_throws DomainError ema_std(RNDS, 3)
    @test_throws DomainError ema_std(RNDS, 5; init_sig=-1.0)

    @test_throws DomainError ema_stats(RNDS, 1)
    @test_throws DomainError ema_stats(RNDS, 2)
    @test_throws DomainError ema_stats(RNDS, 3)
    @test_throws DomainError ema_stats(RNDS, 5; init_sig=-1.0)

    # sig_cumsum input contract tests
    t_test = collect(1.0:20.0)
    @test_throws DomainError sig_cumsum(t_test[1:1], RNDS[1:1], 2, 0.5; chk_inp=true)   # n < 2
    @test_throws DomainError sig_cumsum(t_test[1:5], RNDS, 2, 0.5; chk_inp=true)          # length mismatch
    @test_throws DomainError sig_cumsum(t_test, RNDS, 1, 0.5; chk_inp=true)               # w <= 1
    @test_throws DomainError sig_cumsum(t_test, RNDS, 2, 0.0; chk_inp=true)               # h <= 0
    @test_throws DomainError sig_cumsum(reverse(t_test), RNDS, 2, 0.5; chk_inp=true)      # non-increasing t

    # tic_diff1 / tic_diff2 input contract tests
    @test_throws DomainError tic_diff1(t_test[1:5], RNDS, chk_inp=true)                   # length mismatch
    @test_throws DomainError tic_diff1(reverse(t_test), RNDS, chk_inp=true)                # non-increasing t
    @test_throws DomainError tic_diff2(t_test[1:5], RNDS, chk_inp=true)                   # length mismatch
    @test_throws DomainError tic_diff2(reverse(t_test), RNDS, chk_inp=true)                # non-increasing t

    # entropy_index input contract tests
    @test_throws DomainError entropy_index(RNDS2; n=2)                                     # n <= 2
    @test_throws DomainError entropy_index(RNDS2; λ=0.0)                                   # λ <= 0
    @test_throws DomainError entropy_index(RNDS2; λ=1.1)                                   # λ > 1

    # pow_n input contract tests
    @test_throws DomainError pow_n(2, -1)
    @test_throws DomainError pow_n(2, -1, 5)

    # ewt_mean input contract tests
    ts_test = collect(1.0:20.0)
    xs_test = rand(20)
    @test_throws DomainError ewt_mean(ts_test[1:5], xs_test, 3, 1.0)   # length mismatch
    @test_throws DomainError ewt_mean(ts_test, xs_test, 0, 1.0)         # b <= 0
    @test_throws DomainError ewt_mean(ts_test, xs_test, 20, 1.0)        # b >= n
    @test_throws DomainError ewt_mean(ts_test, xs_test, 3, 0.0)         # lm <= 0
    @test_throws DomainError ewt_mean(ts_test, xs_test, 3, 1.1)         # lm > 1
end

@testset "Finance (tic_diff1/tic_diff2) " begin
    # Test against known polynomial: x(t) = t^2
    # First derivative = 2t, second derivative = 2
    t = collect(1.0:0.5:10.0)
    x = t .^ 2

    d1 = tic_diff1(t, x)
    # Central difference for t^2 should give exactly 2*t at interior points
    for i in eachindex(d1)
        @test d1[i] ≈ 2.0 * t[i+1] rtol = 1e-12
    end

    d2 = tic_diff2(t, x)
    # Second derivative of t^2 = 2
    for i in eachindex(d2)
        @test d2[i] ≈ 2.0 rtol = 1e-10
    end

    # With chk_inp=true, valid inputs should not throw
    d1c = tic_diff1(t, x; chk_inp=true)
    d2c = tic_diff2(t, x; chk_inp=true)
    @test d1c ≈ d1
    @test d2c ≈ d2

    # Test with irregular spacing (central difference gives t_{i+1}+t_{i-1} for t^2)
    t_irreg = [1.0, 1.3, 2.0, 3.5, 4.0]
    x_irreg = t_irreg .^ 2
    d1_irreg = tic_diff1(t_irreg, x_irreg)
    for i in eachindex(d1_irreg)
        @test d1_irreg[i] ≈ t_irreg[i+2] + t_irreg[i] rtol = 1e-12
    end
end

@testset "Finance (sig_cumsum)          " begin
    # A constant series should produce no signals (no deviations from mean).
    t_const = collect(1.0:20.0)
    x_const = fill(5.0, 20)
    tics, sigs = sig_cumsum(t_const, x_const, 3, 0.1)
    @test length(tics) == 0
    @test length(sigs) == 0

    # A step function should produce signals after the step.
    x_step = vcat(fill(0.0, 10), fill(10.0, 10))
    t_step = collect(1.0:20.0)
    tics_s, sigs_s = sig_cumsum(t_step, x_step, 5, 1.0)
    @test length(tics_s) == length(sigs_s)
    @test length(tics_s) > 0
    # All detected tics should be at or after the step
    @test all(tics_s .>= 10.0)

    # With chk_inp=true, valid inputs should not throw
    tics_v, sigs_v = sig_cumsum(t_step, x_step, 5, 1.0; chk_inp=true)
    @test tics_v == tics_s
    @test sigs_v == sigs_s
end

@testset "Finance (ewt_mean)            " begin
    # Constant series: moving average of constant should be constant.
    ts = collect(1.0:50.0)
    xs = fill(3.0, 50)
    wm = ewt_mean(ts, xs, 5, 1.0)
    @test length(wm) == 45
    @test all(wm .≈ 3.0)

    # With decay factor = 1, temporal weighting only
    Random.seed!(42)
    xs2 = rand(50)
    wm1 = ewt_mean(ts, xs2, 5, 1.0)
    @test length(wm1) == 45

    # With decay < 1, result should differ from non-decayed
    wm_decay = ewt_mean(ts, xs2, 5, 0.9)
    @test length(wm_decay) == 45
    @test !(wm1 ≈ wm_decay)

    # Result length should be n - b
    for b in [3, 5, 10]
        wm_b = ewt_mean(ts, xs2, b, 1.0)
        @test length(wm_b) == 50 - b
    end
end

@testset "Finance (PropCheck Test)      " begin

    function prv(v)
        n = length(v)
        idx = Int(floor(rand() * n + 1))
        v1 = reverse(view(v, 1:idx))
        v2 = reverse(view(v, (idx+1):n))
        return (v == reverse(vcat(v2, v1)))
    end

    function permargs(f::Function, w::Vector{T}, tol::T)::Bool where {T<:Real}
        return ((≈)(f(w), f(Random.shuffle(w)); rtol=tol))
    end

    function normalize(w::Vector{T}) where {T<:Real}
        N = length(w)
        nw = copy(w)
        s = zero(T)
        for i in 1:N
            s += nw[i]
        end
        nw ./= s
        return (nw)
    end
        ivec = PropCheck.vector(PropCheck.iconst(N), PropCheck.isample(0.0001:0.001:1.0))
        @test PropCheck.check(v -> permargs(Finance.WWsum, v, TOL2), ivec)
end



