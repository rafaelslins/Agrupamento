
******************************************************************
PROGRAM: 		PPCLUST.SAS
DESCRIPTION:	CLUSTERING BASED ON RANK F TEST OF WANG (2004)
VERSION:		1.0.0
UPDATE:			2018-05
SYSTEM:			SAS WINDOWS 9.4(M4)
MODULES:		BASE - IML - STAT
MACRO:          PPCLUST
NOTE:			ADJUSTED VERSION OF PPCLUST (SORTCENTER)
     			FROM GEORGE FREITAS VON BORRIES

******************************************************************;

libname dados 'C:\Users\Rafael\Desktop\Mestrado\Dissertação\Dados\Alzheimer\data';


%macro ppclust(dataset,obsmin,obsmax,alpha);

   options nosource nosource2 nonotes nostimer nomprint nosymbolgen nomlogic;

   data datanew;
   	 set &dataset;
     idobs = _N_;
   run;

	proc iml;

 * Routine to rank matrix even if it has missing (unbalanced) data;

	start rankmiss(matrix1);
		matrixm = matrix1 = .;
      nmiss = matrixm[,+]; nmiss = nmiss[+,];
      miss = -999999999999;
      matrix1 = choose(matrix1=.,miss,matrix1);
      matrix1 = ranktie(matrix1);
      matrix1 = matrix1 - nmiss;
      matrix1 = choose(matrix1<=0,.,matrix1);
	finish;

* Routine to sort lines (genes) of matrix according to their medians;

	start sortmed(matrix1,matrix2);
   	matrix3 = choose(matrix2=.,max(matrix1)+1000,matrix2);
      medians = t(median(t(matrix3)));
      call sortndx(matrix1,medians,{1});
		matrix2 = matrix2[matrix1,];
      free medians;
   finish;

	start sortcenter(matrix2,colm5,tst);
	   nr = nrow(matrix2);
		nc = ncol(matrix2);
		if tst = 1 then do;
		   m2 = matrix2;
		end;
		else if tst > 1 then do;
			ed = tst-1;
   		m1 = matrix2[1:ed,];
			m2 = matrix2[tst:nr,];
		end;
		nr2 = nrow(m2);
	 	j1 = j(nr2,1) * 0;
		c1 = int(0.40*nr2);
		c2 = int(0.60*nr2);
		*call sort(m2,colm5);
		do i = 1 to nr2;
		  	if i < c1 then j1[i,] = 1;
			if i > c2 then j1[i,] = 1;
		end;
		m2 = m2 || j1;
		nc2 = nc + 1;
		cols = nc2 || nc;
		call sort(m2, cols);
        m2 = m2[,1:nc];
		if tst = 1 then matrix2 = m2;
		else if tst > 1 then matrix2 = m1 // m2;
	finish;

* Routine to calculate sigma4hat - jacknife estimate. Exact estimate
                can be seen in previous versions;

	start sig4hat(x,n,s4hat);
   	mx = sum(x)/n;
      s4hat = (ssq(x-mx)/n)**2;
   finish;

   start sig4jack(x,sigma4);
   	n = nrow(x); ni = n-1;
      if ni > 0 then do;
      	call sig4hat(x,n,s4hat);
         s4jack = n * s4hat;
         do  i=1 to n;
				xni = t(remove(x,i));
            call sig4hat(xni,ni,s4hatni);
            s4jack = s4jack - (n-1)/n * s4hatni;
            free s4hatni;
      	end;
         sigma4 = s4jack;
   	end;
      else if ni = 0 then sigma4 = 0;
   finish;

   start sigma4est(exp1,nrowe1,maxobs,sigma4e);
   	exp1 = colvec(exp1);
      lfac = t(1:nrowe1) @ j(maxobs,1);
      expi1 = exp1 ^= .;
      lfac = lfac # expi1;
      lfac = lfac[loc(lfac),];
      exp1 = exp1[loc(exp1 ^= .),];
      y = j(nrow(exp1),1);
      sigma4e = j(nrowe1,1);
      do i = 1 to nrowe1;
      	temp = choose(lfac=i,y=1,y=0) # exp1;
         temp = temp[loc(temp),];
         call sig4jack(temp,sigma4);
         sigma4e[i,] = sigma4;
      end;
      free exp1 expi1 lfac y temp;
	finish;

* Main Routine for one-way anova based on ranks when we have
  a large number of factors;

   start anovarank(resp,nfac,maxobs,sigma4ea,pv);
   	respi = resp ^= .;
      ng = respi[,+];
      meang = resp[,+]/ng;
      center= resp - repeat(meang,1,maxobs);
      ss = center[,##];
      varg = ss / (ng-1);
      resp = colvec(resp);
      j = j(maxobs,1);
      free respi;
      meant = sum(meang)/nfac;
      difmt = meang - meant;
      mst = ssq(meang - meant)/(nfac - 1);
      mse2= sum(varg/ng) / nfac;
      fr2 = mst/mse2;
      tau2 = sum(sigma4ea/(ng#(ng-1))) * 2 / nfac;
      asyvar = tau2 / mse2**2;
      pv = 1 - probnorm(sqrt(nfac)*(fr2-1)/sqrt(asyvar));
	finish;

   start indtest(exp,st,limitn,alpha,nrowe,ncole,colgr,colts,cols4,g);
	   exp[,colts] = 0; count = 0;
	  	do i=1 to nrowe;
	 		if exp[i,colgr]=g then exp[i,colts] = 1;
		end;
		do i=st to nrowe;
			exp[i,colts] = 1;
			n = sum(exp[,colts]);
	     	expt = exp[loc(exp[,colts]),];
			call anovarank(expt[,1:ncole],n,ncole,expt[,cols4],pv);
			if (pv > alpha) then exp[i,colgr] = g;
			else if (pv <= alpha) then exp[i,colts] = 0;
		end;
		call sort(exp,colgr); * Included;
		do i=1 to nrowe;
		   if exp[i,colgr] <= g then count = count + 1;
		end;
		st = count + 1;
	finish;

	use datanew;
   	read all var{idobs} into factor;
      read all var("&obsmin":"&obsmax") into expression;
      alpha = &alpha;
      testpp = 0; g = 1;
		ncole = ncol(expression); 
		nrowe = nrow(expression);
	  	call rankmiss(expression);
		call sortmed(factor,expression);
		*medians;
		matrix3 = choose(expression=.,max(factor)+1000,expression);
        medians = t(median(t(matrix3)));
		********;
		exp1 = expression;
		call sigma4est(exp1,nrowe,ncole,sigma4e);
		colid = ncole + 1;
		colgr = ncole + 2; 
		colts = ncole + 3;
		cols4 = ncole + 4;
		colm5 = ncole + 5;
		j1 = j(nrowe,1);
		expression = expression || j1 # factor;
		expression = expression || j1 * 9999;
		expression = expression || j1 * 0;
		expression = expression || j1 # sigma4e;
		expression = expression || j1 # factor;
		expression = expression || j1 # medians;
		tstart = 1; tend = nrowe; 
		call sortcenter(expression,colm5,tstart); 
		*print expression;
		call anovarank(expression[,1:ncole],nrowe,ncole,expression[,cols4],pv);	 
		pv2 = pv; 
      if (pv > alpha) then expression[,colgr] = 1;
		if (pv <= alpha) then do;
			L = (tend - tstart + 1)/2;
			tend = tstart + int(L-1);
			do while (tstart <= nrowe);
			   if (pv2 > alpha) then g = g + 1;
			   expression[,colts] = 0;
			   do i = 1 to nrowe;
				   if (tstart <= i) then do;
 						if ( i <= tend) then expression[i,colts] = 1;
					end;
				end;
				n = sum(expression[,colts]); 
			   if (n = 1) then do;
					expression[tstart,colgr] = 0;
					tstart = tend + 1;
					tend = nrowe;	
					testpp = -1;	
				end;
				exp = expression[loc(expression[,colts]),];
				if testpp = 0 then call anovarank(exp[,1:ncole],n,ncole,exp[,cols4],pv);
				pv2 = pv;
				if (pv > alpha) then do;
				   do i = 1 to nrowe;
						if (tstart <= i) then do;
 							if ( i <= tend) then expression[i,colgr] = g;
						end;
					end;
					tstart = tend + 1;
					tend = nrowe;
					st = tstart;  exp = expression; 
					call indtest(exp,st,limitn,alpha,nrowe,ncole,colgr,colts,cols4,g);
					tstart = st;
					expression = exp;
				end;
				else if (pv <= alpha) then do;
					L = (tend - tstart + 1)*0.9;    
					tend = tstart + int(L-1);	
					testpp = 0;
				end;
			end;	
		end;
 
		exp = expression[,colid];
		exp = exp || expression[,colgr];

		create groupclass var{idobs group};
		append from exp;

		free all;
	quit;

   proc sort data=datanew; by idobs;
   proc sort data=groupclass; by idobs; run;

   data groupclass groupclass0;
      merge groupclass datanew;
      by idobs;
      if group = 0 then output groupclass0;
      else if group ^= 0 then output groupclass;
   run;

   proc sort data=groupclass; by group; run;

   data t1 (keep = idobs group &obsmin - &obsmax);
   	set groupclass;
   run;

   data t1 (drop = &obsmin - &obsmax);
   	set t1;
      array num(*) &obsmin - &obsmax;
      length obst $ 8;
      do i=1 to dim(num);
      	expression = num[i];
         call vname(num[i],obst);
         output;
    	end;
   run;

   proc means data=t1 noprint;
   	var expression;
      by group;
      output out=t2 median=median;
   run;

   data t2 (keep=group median);
   	set t2;
   run;

   proc sort data=t2; by median; run;

   data t2 (keep=group corder);
   	set t2;
      corder = _N_;
   run;

   proc sort data=t2; by group; run;

   data groupclass(drop=corder);
   	merge groupclass t2;
      by group;
      group = corder;
	run;
		
	data groupclass; 
  		set groupclass;
  		if idobs=lag(idobs) then delete;
	run; 

	data groupclass;
   	 set groupclass0 groupclass;
    run;

    proc datasets nolist;  
		delete t1 t2 groupclass0;  
	run;

	proc freq data=groupclass noprint;
  		tables group / out=grfr;
	run;
    
    proc sort data=datanew; 
	 	by idobs;
	run;

	proc sort data=groupclass(keep=idobs group);
	    by idobs;

	data datanew;
	  merge datanew groupclass;
	  by idobs;
	run;

	proc sort data=groupclass;
		by group;
	run;

%mend ppclust;			


%ppclust(dados.VCX_dif_bruta,VCX_affected_1,VCX_affected_19,1e-8);
%ppclust(dados.VCX_dif_FC,VCX_affected_1,VCX_affected_19,1e-8);
