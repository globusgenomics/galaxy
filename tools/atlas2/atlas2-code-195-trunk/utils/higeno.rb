
maa=0.99
mab=0.60
mbb=0.10
fst=0.13
pre=10.0

# s is the array of reference and alternative base coverages [ r1, a1, r2, a2, r3, a3...], n is the number of samples
def higeno( s, n )
	normal = Array.new
	n.times do |i|
		ta=s[i*2]
		tb=s[i*2+1]
		normal[i]=Math.lgamma(ta+tb+1)-Math.lgamma(ta+1)-Math.lgamma(tb+1)
	end
	fre =0.0
	old =0.0
	eaa=0.99
	eab=0.6
	ebb=0.10
	naa=1.0
	nab=1.0
	nbb=1.0
	laa=0.0
	lab=0.0
	lbb=0.0
	aap=0.0
	aaq=0.0
	abp=0.0
	abq=0.0
	bbp=0.0
	bbq=0.0
	begin
		laa=Math.log((naa+0.5)/N);	lab=Math.log((nab+0.5)/N);	lbb=Math.log((nbb+0.5)/N);
		aap=Math.log(eaa);	abp=Math.log(eab);	bbp=Math.log(ebb);
		aaq=Math.log(1.0-eaa);	abq=Math.log(1.0-eab);	bbq=Math.log(1.0-ebb);
		naa=nab=nbb=0.0;
		paa=0.0; pab=0.0 ; pbb=0.0; qaa=0.0; qab=0.0 ; qbb = 0.0
		p=s
		n.times do |i|
			ta=p.shift
			tb=p.shift
			nor=normal[i]
			taa=Math.exp(laa+ta*aap+tb*aaq+nor)
			tab=Math.exp(lab+ta*abp+tb*abq+nor)
			tbb=Math.exp(lbb+ta*bbp+tb*bbq+nor)
			ss=1.0/(taa+tab+tbb);	taa*=ss;	tab*=ss;	tbb*=ss;
			paa+=taa*ta;	pab+=tab*ta;	pbb+=tbb*ta;
			qaa+=taa*tb;	qab+=tab*tb;	qbb+=tbb*tb;
			naa+=taa;	nab+=tab;	nbb+=tbb;
		end
		old=fre
		fre=2*naa+nab
		eaa=(paa+maa*pre)/(paa+qaa+pre);
		eab=(pab+mab*pre)/(pab+qab+pre);
		ebb=(pbb+mbb*pre)/(pbb+qbb+pre);
	end while((fre-old).abs>1e-3)

	p=s
	n.times do |i|
		ta=p.shift
		tb=p.shift
		nor=normal[i]
		if(imputation)
			taa=exp(ta*aap+tb*aaq+nor);
			tab=exp(ta*abp+tb*abq+nor);
			tbb=exp(ta*bbp+tb*bbq+nor);
		else
			taa=exp(laa+ta*aap+tb*aaq+nor);
			tab=exp(lab+ta*abp+tb*abq+nor);
			tbb=exp(lbb+ta*bbp+tb*bbq+nor);
		end
		ss=1.0/(taa+tab+tbb);	taa*=ss;	tab*=ss;	tbb*=ss;
		if(taa>tab&&taa>tbb)	
			print "0/0:#{tb}:#{ta}:#{ta+tb}:#{(-10*Math.log10([1-taa, 1e-3].max)).round}\t";
		elsif(tab>taa&&tab>tbb)	
			print "0/1:#{tb}:#{ta}:#{ta+tb}:#{(-10*log10([1-tab, 1e-3].max)).round}\t";
		else	
			print "1/1:#{tb}:#{ta}:#{ta+tb}:#{(-10*log10([1-tbb,	1e-3].max)).round}\t";
		end
	end
	put "\n"
end
