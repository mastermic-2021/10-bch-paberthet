/* x is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
f = ffgen(q,'a)
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);

/*On va utiliser le décodage par syndrome*/
\\commençons par traduire le chiffré en forme polynomiale pour obtenir y(x)

y = int2fqx(m,f);

\\on calcule ensuite le syndrome :

s(y,alpha,l) = sum(i=0,2*t-1, subst(y, 'x, alpha^(i+l))*'x^i);

\\on emploie ensuite les approximants de Padé pour trouver le message d'origine : 

correction(n) = {
for(k = 0, q, alpha = ffprimroot(f);
for(j = 0, q-1, err = List(); pade = bestapprPade(Mod(s(y,alpha,j),x^(2*t))); R = numerator(pade); E = denominator(pade);
for(i = 0, q-2, P = subst(E,'x,alpha^(-i)); if(P==0, val = subst(R/deriv(E) * (x^(j-1)), 'x , alpha^(-i));list = fqx2int(val,f); listput(err, list)));
if(#err >10, print(Strchr(Vec(err)));return);
);
);
};

correction();
