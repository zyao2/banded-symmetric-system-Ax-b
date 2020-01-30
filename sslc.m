function  x=sslc(n,n2,A,beg,nd)
%!__________________________________________________________________
%!
%!    solving linear symmetrical system Ax=b with several diagonales matrix
%!__________________________________________________________________
%!    Input DATA 
%!         n   - number of unknowns
%!         n2  = 2*n
%!         nd  - number of diagonales
%!         A   - matrix with diagonales  (n*nd )     
%!         beg - right part of the system (first n)
%!         x   - solution of the system (first n)
%!         w   - working array(n2,nd)
%!         y   - working array(n2)
%!__________________________________________________________________
%!        Diagonals outside matrix should be zeroed!!!
y = zeros(n2,1);
w = zeros(n2* nd,1);

nd1=nd-1;
EPD=0.00001;
nnd=n+nd;
x=zeros(n,1);
  
if(n == 1) 
	q=A(1);
	if(fabs(q) < EPD) 
        q=EPD;
    end
	x(1)=beg(1)/q;
	return;
end
if(n == 2)
    A11=A(1);
	A21=A(2);
	A12=A(3);
	q=A11*A21 - A12*A21;
	if(abs(q) < EPD)
        q=EPD;
    end
	B1=beg(1);
	B2=beg(2);
	x(1)=(B1*A21-B2*A12)/q;
	x(2)=(B2*A11-B1*A12)/q;
	return ;
end
  
for i=1:nd1
    x(nnd-i)=0.;
	y(i)=0.;
	for j=1:nd
        w(i+(j-1)*n2)=0.;
    end
end
for j=2:nd
    im=j-1;
    for i=1:im
        A((j-1)*n+n-i+1)=0;
    end
end

for i=1:n
    iw=i+nd1;
	y(iw)=beg(i);
	div=1.;
	j=1;
    s=A(i+n*(j-1));
	km=nd-j;
	if (j ~= nd ) 
		for k=1:km
            nw=nd-k;
            lw=i+j+k-2;
            nw1=nw-j+1;
			ii=lw+nw*n2;
			s=s-w(ii)*w(lw+nw1*n2)*w(lw);
			y(iw)=y(iw)-w(ii)*y(lw);
        end  
    end
	w(iw)=s/div; 
    if(abs(w(iw)) < EPD)  
        w(iw)=EPD;
    end

    div=w(iw);

	for j=2:nd
		s=A(i+n*(j-1));
		km=nd-j;
		if (j ~= nd )
            for k=1:km
                nw=nd-k;
				lw=i+j+k-2;
				nw1=nw-j+1;
				ii=lw+nw*n2;
				s=s-w(ii)*w(lw+nw1*n2)*w(lw);
            end
        end
		w(iw+(j-1)*n2)=s/div; 
    end
end
for k=1:n;
    i=n+1-k;
	l=n+nd-k;
	x(i)=y(l)/w(l);
	for j=1:nd1
        x(i)=x(i)-w(l+j*n2)*x(i+j);
    end
end

