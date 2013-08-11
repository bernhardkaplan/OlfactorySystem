NEURON {
	SUFFIX odorinput
	NONSPECIFIC_CURRENT i
	RANGE i, e, gor, tstart,or
}

PARAMETER {
	or = 0.0 <0,1>
	gor = 0.001 (siemens/cm2) <0, 1e9>
	e = 0 (millivolt)
	tstart = 0  (ms) :will be overwritten by cell constructor parameters
	tau = 20 (ms) 
	tstop = 500 (ms) :will be overwritten by cell constructor parameters
	: good value
	:tstop = tstart + 25 * tau(ms) :will be overwritten by cell constructor parameters
	: these values give t_rise = 87 ms, t_stim =372 ms
}

ASSIGNED {
	i (milliamp/cm2)
	v (millivolt)
}

BREAKPOINT {
	i = 0
		
	:if (t>=tstart && t<=tstop){ 	: stimulus going in
	:	i = or * gor * (v - e)
	:	i = (sin(3.141593*t/278))^2 * or * gor * (v - e)
	:if((t-tstart)<200){			: if we are in the first 200ms of the stimulus, it is regulated by a sigmoid function
	:	i = (1/(1+exp(-1/20*((t-(tstart+100)))))) * or * gor * (v - e)
			: fact = ((t-tstart)/200)
	:	}
	:}
	

	if (t>=tstart && t<=tstop){ 	: stimulus increase
		if (t > tstart + 30 * tau) {
			i = or * gor * (v - e)
		}
		else {
			i = (1 / (1 + exp(-(t - 12*tau) / tau))) * or * gor * (v - e)
		}
	}
	if (t > tstop){ : stimulus decrease
        i = (1 / (1 + exp(-(tstop + 10 * tau - t) / tau))) * or * gor * (v - e)
	}
	: else i = 0
}

