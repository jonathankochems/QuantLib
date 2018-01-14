import QuantLib as ql # version 1.5

s0 = 100.
v0 = 0.04
alpha = 1.0
b = 0.04
sigma = 2.0
rho = -0.7
Tdays = 30
K = 100.

daynow = 0
day_count = ql.Actual365Fixed()
time_zero = ql.Date(27,7,2017)
calculation_date = time_zero + daynow
ql.Settings.instance().evaluationDate = calculation_date

spot_handle = ql.QuoteHandle( ql.SimpleQuote(s0))
dividend_rate =  0.0
risk_free_rate = 0.0
flat_ts = ql.YieldTermStructureHandle( ql.FlatForward(calculation_date, risk_free_rate, day_count) )
dividend_yield = ql.YieldTermStructureHandle(ql.FlatForward(calculation_date, dividend_rate, day_count))
heston_process = ql.HestonProcess(flat_ts,dividend_yield,spot_handle,v0,alpha,b,sigma,rho)


# option data
maturity_date = time_zero + Tdays
strike_price = K
payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike_price)
exercise = ql.EuropeanExercise(maturity_date)
european_option = ql.VanillaOption(payoff, exercise)

engine = ql.AnalyticHestonEngine(ql.HestonModel(heston_process),128)
european_option.setPricingEngine(engine)
heston_price = european_option.NPV()
print "The theoretical price is ", heston_price
