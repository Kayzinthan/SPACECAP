SPACECAP <- function()
{
    # Loads tcktk package, TeachingDemos and the tck package Tktable
    require(tcltk) 
    library(TeachingDemos) 
    tclRequire("Tktable") 
  
    # Create a top level window from a tkwidget
  	tt <- tktoplevel()
  
  	# Define Global variables 
  	locidso <- NULL
  	locso <- NULL
  	grid2500 <- NULL
  	Global_iterNum <- 2500
  	Global_burn <- 500
  	Global_skip <- 4
  	Global_nz <- 250
  	Global_acSize <- 0.336
    statusText <- tclVar(" ")
    shouldIStop <- FALSE
  
    ## Action function to read animal capture details file
  	readFile1 <- function ()
  	{
  		file <- tclvalue(tkgetOpenFile())
   		if (!length(file)) return()
  		locidso <<- read.csv(file)
  		tkconfigure(capt, state="normal")
		  tkdelete(capt, "0", "insert")
  		tkinsert(capt, "0", file)
  		tkconfigure(capt, state="disable")
  		tkfocus(tt)
  	}
  
    ## Action function to read trap deployment details file
  	readFile2 <- function ()
  	{
  		file <- tclvalue(tkgetOpenFile())
  		if (!length(file)) return()
  		locso <<- read.csv(file)
  		tkconfigure(trap, state="normal")
		  tkdelete(trap, "0", "insert")
  		tkinsert(trap, "0", file)
  		tkconfigure(trap, state="disable")
  		tkfocus(tt)
  	}
  
    ## Action function to read potential activity centers file
  	readFile3 <- function ()
  	{
  		file <- tclvalue(tkgetOpenFile())
  		if (!length(file)) return()
  		grid2500 <<- read.csv(file)
  		tkconfigure(grids, state="normal")
		  tkdelete(grids, "0", "insert")
  		tkinsert(grids, "0", file)
  		tkconfigure(grids, state="disable")
  		tkfocus(tt)
  	}

    ## Define the frame structure for SPACECAP
    one <- tkframe(tt, width = 800, height = 250, bg= "light blue",relief="groove",borderwidth=2)
    tkpack(one, side = "top")
    two <- tkframe(tt, width = 800, height = 300, bg= "white",relief="groove",borderwidth=2)
   	tkpack(two, side = "top")
    three <- tkframe(tt, width = 800, height = 40, bg= "white",relief="groove",borderwidth=2)
    tkpack(three, side = "bottom")

    # Sub-Frames
  	fileFrame <- tkframe(one, width = 300, height = 250, bg= "light blue",relief="groove",borderwidth=2)
  	tkpack(fileFrame, side = "left")
    choiceFrame <-tkframe(one, width = 250, height = 250, bg= "light blue",relief="groove",borderwidth=2)
    tkpack(choiceFrame, side = "left")
    valueFrame <-tkframe(one, width = 250, height = 250, bg= "light blue",relief="groove",borderwidth=2)
    tkpack(valueFrame, side = "right")
    
    # Define status window
    statusWin <- tktext(three,height=5)
    scr <- tkscrollbar(three, command=function(...) tkyview(statusWin,...))
    tkconfigure(statusWin, yscrollcommand=function(...) tkset(scr,...))
    tkpack(scr, side="right", fill="y")
    tkpack(statusWin, fill="both", expand=TRUE)
     
    ## Function for reading input data and creating 3D data structures from it
    readData <- function()
  	{
    		nID = max(locidso[,2])
    		nSO = (dim(locso)[2]) - 3
    		nLOC = dim(locso)[1]
        ### Function for getting capture values in a 3 dimensional ID x SO x LOC format
    		makeData3d = function()
    		{
          data3d = structure(rep(0, times=nID*nSO*nLOC), .Dim=c(nID, nSO, nLOC))
          len = length(locidso[,1])
          for (i in 1:len)
          {
            loc = locidso[i,1]
            id = locidso[i,2]
            so = locidso[i,3]
            data3d[id,so,loc] = 1
          }
    		  data3d
    		}
        ### Function for getting the deployment values in a LOC x SO format
    		makeMask3d = function()
    		{
    		  ## In the deployment file, the SO start from column 4
      		posFirstSO = 4
      		posLastSO = nSO+3
          mask3d = as.matrix(locso[,posFirstSO:posLastSO])

      		mask3d
    		}
        tiger3dData <- structure(list(makeData3d(), makeMask3d()), .Names = c("data3d", "mask3d"))
    		grid2500 <- structure(c(grid2500[,1], grid2500[,2], grid2500[,3]), .Dim = c(length(grid2500[,1]), 3), .Dimnames = list(c(1:length(grid2500[,1])), c("X_Coord", "Y_Coord", "HABITAT")))
    		ctLocs <- structure(list(grid = structure(list(x = locso[,2], y = locso[,3]), .Names = c("x", "y"), row.names = c(1:nLOC), class = "data.frame")), .Names = c("grid"))

        list(tiger3dData=tiger3dData, grid2500=grid2500, ctLocs=ctLocs)
    } ## End of readData function
    	
    ####
    ## The main script for running the Sapcecap analysis
    ####
    SCRd.fn <- function(ni=52000,burn=2000,skip=50,bsigma=1,Mb=0,nz=450,dexp=2)
    {
        # Feb 23 -- andy tinkered with activity center updating
        # Feb 24 -- changed behavioral response to regression variable
        # Sept 26 -- work on including Mb.
        
        # ni = number of iterations, total
        # burn = number to discard
        # skip = thin rate (i.e., keep one in skip iterations)
        # bsigma = 0 fits non-spatial model
        # GRID is the discrete grid of the study area
        # nz = number of "all zero" encounter histories to add
        
        
        # Required data objects are the encounter array Y (bindary observations)
        # which is nind x T x ntraps
        # MASK which is ntraps x reps (it gets transposed though)
        # traplocs = coordinates of trap locations

#################### GUI related START######################

      	cat("Starting Analysis...", fill=TRUE)
      	
        statusText <<- "\nStarting analysis\n"
        tkinsert(statusWin, "end", statusText)
	
#################### GUI related END #######################		
		
      	# Read data from input files 
      	hold = readData()
    
      	tiger3dData = hold$tiger3dData
      	grid2500 = hold$grid2500
      	ctLocs = hold$ctLocs
      
      	Y<-tiger3dData$data3d          # nind x nT x ntraps
      	MASK<-t(tiger3dData$mask3d)    # want this to be rep x traps
      	traplocs<-as.matrix(ctLocs$grid)
      	GRID=grid2500
		
		
################### Output related START ####################
		gridsub=subset(grid2500, grid2500[,3]>0)
		nhrc=sum(grid2500[,3])
		

################### Output related END ######################
		
        ###
        ###
        # what are the dimensions of the problem
        ###
        nind<-dim(Y)[1]
        nT<-dim(Y)[2]
        M<-nind+nz   # total size of data set
        ntraps<-nrow(traplocs)
  
  
        ###
        ###
        ## following lines scale coordinate system. This is kind of arbitrary.
        ###
        ###
        mgx<-min(traplocs[,1])
        mgy<-min(traplocs[,2])
        traplocs[,1]<-(traplocs[,1]-min(traplocs[,1]))/5000
        traplocs[,2]<-(traplocs[,2]-min(traplocs[,2]))/5000
        G<-GRID[,1:2]
        G[,1]<-(G[,1]-mgx)/5000
        G[,2]<-(G[,2]-mgy)/5000
        goodbad<-GRID[,3]
        G<-G[goodbad==1,]
        nG<-nrow(G)
  
  
        ###
        ### Data processing -- this block of code determines a neighborhood
        ### for every pixel. That information is used in the MCMC updating
        ### of activity centers. Some tuning would be required
        NN<-matrix(NA,nrow=nG,ncol=400)
        # RAD is a number related to the MCMC for the activity centers.  It should be
        # big enough so that every grid point has "several" neighbors
        # for 10k grid I used RAD=1
        RAD<-1
        for(i in 1:nG){
        od<- sqrt( (G[i,1]-G[,1])^2  +  (G[i,2]-G[,2])^2  )
        od<- (1:length(od))[od<RAD]
        NN[i,1:length(od)]<-od
        }
        numnn<-apply(!is.na(NN),1,sum)
        NN<-NN[,1:max(numnn)]
        #print(min(numnn))
        
        
        ###
        ### does the data augmentation
        ###
        Yaug<-array(0,dim=c(nind+nz,nT,ntraps))
        for(j in 1:nind){
        Yaug[j,1:nT,1:ntraps]<-Y[j,1:nT,1:ntraps]
        }
        M<-nind+nz
  
  
        `e2dist` <- function (x, y)
       	{
        		i <- sort(rep(1:nrow(y), nrow(x)))
        		dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
        		matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
     	  }

        
        ###
        ###
        ## this block of code picks starting coordinates for each individual
        ###
        ###
        
        centers1<-rep(NA,nind)
        for(i in 1:nind){
        tt<-Yaug[i,,]
        tt<-col(tt)[tt==1]   # which traps was he captured in 
        xxx<-traplocs[tt,]  ## coordinates of those traps
        av.coord<-apply(matrix(xxx,ncol=2),2,mean) ## takes average coordinate
        dvec<-as.vector(e2dist(matrix(av.coord,ncol=2),G))  # finds closest grid pt
        centers1[i]<-(1:length(dvec))[dvec==min(dvec)][1]   # that is initial loc
        }
        # uncaptured guys need centers too.....
        centers2<-sample(1:nG,M-nind,replace=TRUE)
        centers<-c(centers1,centers2)
        S<-G[centers,]   # initial locations for all M individuals
  
  
        ###
        # create "Data" vector but with trap mask information
        ###
        msk<-MASK
        msk2<-array(NA,c(nind+nz,nT,ntraps))
        for(i in 1:(nind+nz)){
        msk2[i,1:nT,1:ntraps]<-msk[1:nT,1:ntraps]
        }
        msk2<-as.vector(msk2)
  
    
        ###
        # create covariate of previous capture
        ###
        prevcap<-array(0,c(nind+nz,nT,ntraps))
        for(i in 1:(nind)){
          for(j in 1:ntraps){
          tmp<-Yaug[i,1:nT,j]
            if(any(tmp==1)){
             fst<- min( (1:nT)[tmp==1] )
             if(fst<nT)
              prevcap[i,(fst+1):nT,j]<-1 
            }
          }
        }
        prevcap<-as.vector(prevcap)


        ##
        ## vectorize all the data objects
        ##
        arr.trues <- array(TRUE, c(nind+nz,nT,ntraps))
        idx<-which(arr.trues, arr.ind = TRUE)
        y<-as.vector(Yaug)
        y<-y[msk2==1]
        prevcap<-prevcap[msk2==1]   #### + 1   ### add 1 if want 1/2 otherwise dummy
        indid<-idx[msk2==1,1]
        repid<-idx[msk2==1,2]
        trapid<-idx[msk2==1,3]
        
  
        ###
        ###
        ## starting values of parameters and other stuff needed, including utility funcs
        ###
        ###
        
        clogloginvpart1<-function(lp){
        # goes from linear predictor, lp, to log(1-p)
         -1*exp(lp)
        }
        clogloginvpart2<-function(part1result){
        # goes from  log(1-pi) to pi
        1-exp(part1result)
        }
  
  
        trapgridbig<-traplocs[trapid,]   # stretches out the trap coord matrix 
        y1<-y==1
        c1<- (S[indid,1]-trapgridbig[,1])^2
        c2<- (S[indid,2]-trapgridbig[,2])^2
  
        sigma<- 5*sqrt(0.30/2)
        loglam0<-log(.018)
        beta<-0
        lam0<-exp(loglam0)
        psi<-.5  # not a good starting values
        delta<- 0 
        z<-c(rep(1,nind),rbinom(nz,1,psi))
  
################ Output related START ###################
        out<-matrix(NA,nrow=(ni-burn)/skip,ncol=5)
        dimnames(out)<-list(NULL,c("sigma","lam0","beta","psi","Nsuper"))
        zout<-matrix(NA,nrow=(ni-burn)/skip,ncol=M)
        Sout<-matrix(NA,nrow=(ni-burn)/skip,ncol=M)
#		indlocs<-matrix(NA,nrow=(ni-burn)/skip,ncol=M)
#		pixcount<-c(rep(0, length(grid2500[ ,2])))
#		pixdensity<-c(rep(0, length(pixcount)))

		
        ################ Utility functions for generating spatial plots #################
		
        ####### Setting the scale for spatial plots #######
#		imagescale <-function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks", "ranges")){
			
			# sort out the location
#			n <- length(col)
#			usr <- par("usr")
#			mx <- mean(usr[1:2]); my <- mean(usr[3:4])
#			dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
#			if (missing(x))
#			x <- mx + 1.05*dx/2	# default x to right of image
#			else if (is.list(x)) {
#				if (length(x$x) == 2) 
#				size <- c(diff(x$x), -diff(x$y)/n)
#				y <- x$y[1]
#				x <- x$x[1]
#			} else x <- x[1]
#			if (is.null(size))
#			if (is.null(y)) {
#				size <- 0.618*dy/n	# default size, golden ratio
#				y <- my + 0.618*dy/2	# default y to give centred scale
#			} else size <- (y-my)*2/n
#			if (length(size)==1)
#			size <- rep(size, 2)	# default square boxes
#			if (is.null(y))
#			y <- my + n*size[2]/2
			
			# draw the image scale
#			i <- seq(along = col)
#			rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2], 
#				 col = rev(col), xpd = TRUE)
			
			# sort out the labels
#			rng <- range(z, na.rm = TRUE)
#			bks <- seq(from = rng[2], to = rng[1], length = n + 1)
#			bks <- formatC(bks, format="f", digits=digits)
#			labels <- match.arg(labels)
#			if (labels == "breaks")
#			ypts <- y - c(0, i) * size[2]
#			else {
#				bks <- paste(bks[-1], bks[-(n+1)], sep = " - ")
#				ypts <- y - (i - 0.5) * size[2]
#			}
#			text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj =
#				 ifelse(size[1]>0, 0, 1), xpd = TRUE) 
#			}
		
		
		####### the spatial plot utility function #######
#		spatialplot<-function(x,y){
#			nc<-as.numeric(cut(y,20))
#			plot(x,pch=" ")
#			points(x,pch=20,col=topo.colors(20)[nc],cex=2)
#			imagescale(y,col=topo.colors(20))
#			}
		
######### Output related END ###############		
######### GUI related START ################		
        
        if(shouldIStop==TRUE) {
            shouldIStop <<- FALSE; 
            statusText <<- "Analysis stopped\n"
            tkinsert(statusWin, "end", statusText)
            stop(call.=FALSE, "Analysis stopped") }
  
      #Progress Bar
	    pb <- tkProgressBar(title = "SPACECAP Progress Bar", min = 0, max = ni, width = 600)

      statusText <<- "Burn-in in progress\n"
      tkinsert(statusWin, "end", statusText)
######### GUI related END ##################
      m<-1
      for(i in 1:ni)
      {
######### GUI related START ##################
          cat("iter: ",i,fill=TRUE)
          Sys.sleep(0.1)
          if(shouldIStop==TRUE) {
            shouldIStop <<- FALSE; 
            statusText <<- "Analysis stopped\n"
            tkinsert(statusWin, "end", statusText)
            stop(call.=FALSE, "Analysis stopped") }
          setTkProgressBar(pb, i, label=paste(round((i*100)/ni),"% completed"))
######### GUI related END ##################            
          ########################
          ########################
          ########################
          ########################
          # PART 1 OF THE MCMC ALGORITHM UPDATES THE REGRESSION PARAMETERS. FOR THIS MODEL
          # THE REGRESSION PARAMETERS ARE (1) INTERCEPT (2) EFFECT OF PREVIOUS CAPTURE
          # (BEHAVIORAL RESPONSE) (3) THE SPATIAL PARAMETER "sigma"
          ### Updating parameters here should only involve guys with z = 1 (i.e., members of the population)
  
          lp<-  loglam0 + Mb*beta*prevcap - ((12.5*bsigma)/(sigma^2))*((c1+c2)^(dexp*0.5))
		  lik1<- log( expm1(exp(lp)))
          lik2<- clogloginvpart1(lp)
          llvector<- lik2
          llvector[y1]<- llvector[y1]+ lik1[y1]
  
          ### generate candidate values
          loglam0c<-rnorm(1,loglam0,.015)
          sigmac<-exp(rnorm(1,log(sigma),.15))
  
          lpc<-  loglam0c + Mb*beta*prevcap - ((12.5*bsigma)/(sigmac^2))*((c1+c2)^(dexp*0.5))
          lik1c<- log(expm1(exp(lpc)))
          lik2c<- clogloginvpart1(lpc)
          llvector.new<- lik2c
          llvector.new[y1]<- llvector.new[y1]+ lik1c[y1]
         
          if(runif(1)<exp(sum( rowsum(llvector.new-llvector,indid)[z==1])))
          {
             loglam0<-loglam0c
             lam0<-exp(loglam0)
             sigma<-sigmac
             llvector<-llvector.new
          }
  
  
          ### Sept 26 2009 added block of code below
          ### to deal with model Mb.
          betac<- rnorm(1,beta,.05)
          lpc<-  loglam0 + Mb*betac*prevcap - ((12.5*bsigma)/(sigma^2))*((c1+c2)^(dexp*0.5))
          #lik1c<- log( (1/exp(-exp(lpc)))  -1 )
          lik1c<- log( expm1(exp(lpc)))
          
          lik2c<- clogloginvpart1(lpc)
          llvector.new<- lik2c
          llvector.new[y1]<- llvector.new[y1]+ lik1c[y1]
          if(runif(1)<exp(sum( rowsum(llvector.new-llvector,indid)[z==1])))
          {
             beta<- betac
             llvector<-llvector.new
          }
          ##############################
          ##############################
          
  
          ########################
          ########################
          ########################
          ########################
          # PART 2 OF THE MCMC ALGORITHM UPDATES THE DATA AUGMENTATION PARAMETERS
          # THIS INCLUDES THE LATENT "z" VARIABLES AS WELL AS THE CRITICAL 
          # PARAMETER "psi"
          ########################
          ########################
          
          # This is the data augmentation part. This code updates each z[i]
          # for all i=1,2,...M individuals. z[i]=1 is a "real" animal that lives
          # in S, whereas z[i]=0 are excess zeros in the data set
          # this is an application of Bayes rule to get Pr(z=1| y[i,,]=0)
          
          
          probz<- exp( rowsum(llvector[indid>nind],indid[indid>nind]) ) # only for nind+1...M
          probz<- (probz*psi )/(probz*psi + (1-psi))
          z[(nind+1):M]<-rbinom(M-nind,1,probz)
          psi<-rbeta(1,1+sum(z),1+M-sum(z))
  
  
          ########################
          ######################## 
          ########################
          ## PART III OF THE ALGORITHM -- THIS BLOCK OF CODE UPDATES THE
          ## ACTIVITY CENTERS
          ###
          # in practice only have to do this for guys with z = 1
          # from the data augmentation.......
          # guys with z = 0 should be drawn uniformly
          # and their acceptance probability should be set to 1.0
          
          newcenters<-trunc(runif(M,0,numnn[centers]))+1  # this might be a dumb way to generate random integers.
          newcenters<- NN[cbind(centers,newcenters)]
          
          # these probabilities are needed in the Metrpolis "acceptance probability" calculation
          qnew<- 1/numnn[centers] 
          qold<- 1/numnn[newcenters]
          
          ## This lines redraws the z=0 ones to be uniform
          ### Can use same proposal as above, no problem 
          Sc<-G[newcenters,]
          c1c<- (Sc[indid,1]-trapgridbig[,1])^2
          c2c<- (Sc[indid,2]-trapgridbig[,2])^2
          
          # note sometimes xx can be 0 which evaluates log(xx/(1-xx)) to -Inf
          # I am not sure if this is really a problem
          xx<-1-exp(-(exp(loglam0+ Mb*beta*prevcap[y1]) )* exp(-((12.5*bsigma)/(sigma^2))*((c1c[y1]+c2c[y1])^(dexp*0.5))))
          xx<- log( xx/(1-xx)  )
          zz<-    -exp(loglam0+Mb*beta*prevcap)*exp(-((12.5*bsigma)/(sigma^2))*((c1c+c2c)^(dexp*0.5)))
          llvector.tmp<- zz
          llvector.tmp[y1]<-llvector.tmp[y1]+xx
          likdiff<-rowsum(llvector.tmp-llvector,indid)
          likdiff[z==0]<-0   # this line was in wrong place if using local proposal as above
          
          # this lines sets acceptance prob to 1 for z=0 guys. Note: some calcs in prev
          # lines probably dont have to be done .... but may be more costly to not do them
          likdiff<-likdiff + log(qold/qnew)
          accept<- runif(M)<exp(likdiff)
          #cat("accept rate: ",mean(accept),fill=TRUE)
          S[accept,]<-Sc[accept,]
          centers[accept]<-newcenters[accept]
          c1<- (S[indid,1]-trapgridbig[,1])^2
          c2<- (S[indid,2]-trapgridbig[,2])^2
          
          ### should update llvector right here.................
          
          ########################
          ########################
          ########################
          
          ##
          ## save output if iteration is past burn-in and using a thin-rate of "skip"
          ##
          
          ## PS : Adding an additional condition to check if skip==1
          if( (i>burn) & ((i%%skip == 1) | (skip ==1)) )
          {
            zout[m,]<-z
            Sout[m,]<- centers
            out[m,]<-c(sigma,lam0,beta,psi,sum(z))
            
             if(Mb == 0)       #Rashmi
             {
             #print ("Line 577: Out matrix contents")
             cat ("sigma", "   lam0", "    psi", "    Nsuper", "\n")
             cat (round(out[m,1], 6), round(out[m,2], 6), round(out[m,4], 6), round(out[m,5], 6), "\n")
             }
             else
             {
             #print ("Line 582: All out matrix contents")
             print (out[m,])
             }
            
            #print(out[m,])
########## GUI related START ############

            if(Mb == 0)     #Rashmi
            {
            statusText <<- paste("sigma\tlam0\tpsi\tNsuper\n",round(out[m,1],2), "\t",round(out[m,2],2), "\t",round(out[m,4],2), "\t",round(out[m,5],2),"\n") 
            }
            else
            {
            statusText <<- paste("sigma\tlam0\tbeta\tpsi\tNsuper\n",round(out[m,1],2), "\t",round(out[m,2],2),"\t", round(out[m,3],2), "\t",round(out[m,4],2),"\t",round(out[m,5],2),"\n")
            }
            
            tkinsert(statusWin, "end", statusText)
########## GUI related END #############
            m<-m+1
          }

      } # End of the iterations loop
		
########## Output related START ###########
      ## Time stamp
      ts <- format(Sys.time(), "%H:%M:%S")
      ts <- unlist(strsplit(ts, ":"))
      ts <- paste(ts[1], ts[2], ts[3], sep="")
      
      #list(out=out,G=G,traplocs=traplocs,Sout=Sout,zout=zout)
      folderName <- paste("output_", ts, sep="")
      dir.create(path=folderName)
      fname <- paste(folderName, "/param_val_", ts,".csv", sep="")
      write.csv(file=fname, out)
      cat("Analysis Complete. Parameter estimates written to ", getwd(), "/", fname, sep="", fill=TRUE)
      
      ## Calculations to obtain pixel-specific densities
#	  indlocs<-Sout*zout
#	  freqonpixel<-as.data.frame(table(indlocs))
	  
#		for (i in 1:length(freqonpixel[ ,1])){
#			if(freqonpixel[i,1]!=0){
#			pixcount[freqonpixel[i,1]]<-freqonpixel[i,2]}
#		}
#		pixdensity<-pixcount/m
      ### logpixdensity<-log(pixdensity)
########### Output related END ############	 
	  
	  ## Addition on July 06 2011 to generate the indicators_val and centers_val output files
#     nameoffile1 = paste(folderName,"/centers_val_",ts,".csv", sep="")
#    	write.csv(Sout, file =nameoffile1)
#	  nameoffile2 = paste(folderName,"/indicators_val_",ts,".csv", sep="")
#	    write.csv(zout, file =nameoffile2)
		
########## Addition on October 21, 2011 to generate pixel densities ########
	  niter=nrow(Sout)
		
# Define a frequency and density matrix
	  freqmat=matrix(data=NA,nrow=nhrc,ncol=2)
	  densitymat=matrix(data=NA, nrow=nhrc, ncol=2)
# Define locations where real animals have been captured
	  indlocs = Sout*zout
# Obtain a count of trap-specific encounters and create pixel-specific densities
		for (i in 1:nhrc){
			freqmat[i,1]=i
			densitymat[i,1]=i
			freqmat[i,2]=length(which(indlocs==i))
			densitymat[i,2]=freqmat[i,2]/niter
	}
#Combining X & Y coordinates to densitymat file
		cdensitymat<-cbind(gridsub[,1:3],densitymat[,2])
		        colnames(cdensitymat)[4]<-'Pixel Density'
   
#Generate csv file for pixel densities
		nameoffile3 = paste(folderName,"/pixeldensities_val_",ts,".csv", sep="")
	    write.csv(cdensitymat, file =nameoffile3)

		
#	  nameoffile4 = paste(folderName,"/freqonpixel_",ts,".csv", sep="")
#	    write.csv(freqonpixel, file =nameoffile4)
#	  nameoffile5 = paste(folderName,"/indlocs_",ts,".csv", sep="")
#	    write.csv(indlocs, file =nameoffile5)
#	  nameoffile6 = paste(folderName,"/pix_densities_", ts, ".csv", sep="")
#		write.csv(pixdensity, file=nameoffile6)
      ##
############ GUI related START #############      
      statusText <<- paste("Analysis Complete\nParameter estimates written to ", getwd(), "/", fname, "\n", sep="") 
      tkinsert(statusWin, "end", statusText)
      close(pb)
############ GUI related END ################   
############ Output related START ##############
      ## Calculating summary statistics of paramteres
      derivedDen <- (out[,5]/(nG*as.numeric(Global_acSize)))*100
      p1 <- 1-exp(-1*out[,2])
      p2 <- 1-exp(-1*out[,3])
      newout <- cbind(out, derivedDen, p1, p2)
      resTable = matrix(NA, nrow=8, ncol=4)
      resTable[,1] = round(apply(newout, 2, mean), digits=4)
   		resTable[,2] = round(apply(newout, 2, sd), digits=4)
       for (i in 1:8){
        resTable[i,3:4] = round(emp.hpd(newout[,i], conf=0.95), digits=4) }

    	rownames=c("_","sigma", "lam0", "beta", "psi", "Nsuper", "Density", "p1", "p2")
      dim(rownames)<-c(9,1)
      colnames=c("Posterior_Mean", "Posterior_SD", "95%_Lower_HPD_Level", "95%_Upper_HPD_Level")
      dim(colnames)<-c(1,4)
      resTable=cbind(rownames,rbind(colnames,resTable))

      if(Mb == 0)   #Rashmi: file has been renamed filename and filename has been renamed as nameoffile to avoid partial argument matching warnings
      {
         for ( i in 1:2 ) {
        nameoffile = paste("./", folderName, "/density_", rownames[i+1,1], "_", ts, ".jpeg", sep="")
        jpeg(filename=nameoffile)
        plot(density(out[,i]), main=rownames[i+1,1])
        dev.off()      
        }
        
        for ( i in 4:5 ) {
        nameoffile = paste("./", folderName, "/density_", rownames[i+1,1], "_", ts, ".jpeg", sep="")
        jpeg(filename=nameoffile)
        plot(density(out[,i]), main=rownames[i+1,1])
        dev.off()      
        }
      }
      else
      {
         for ( i in 1:5 ) {
        nameoffile = paste("./", folderName, "/density_", rownames[i+1,1], "_", ts, ".jpeg", sep="")
        jpeg(filename=nameoffile)
        plot(density(out[,i]), main=rownames[i+1,1])
        dev.off()      
        }
      }

      
      
	  ##### Spatial plot #####
#	  nameoffile = paste("./", folderName, "/SpatialPlot_", ts, ".jpeg", sep="")
#		jpeg(file=nameoffile)
#		spatialplot(gridforplot, pixdensity)
#		dev.off()
############ Output related END ##############	
############ GUI related START ###############		

      if(Mb == 0)   #Rashmi
      {
      statusText <<- paste("Density plots of parameters sigma, lam0, psi, Nsuper saved in jpg format to ", getwd(), "/", folderName, "\n", sep="")
      }
      else
      {
      statusText <<- paste("Density plots of parameters sigma, lam0, beta, psi, Nsuper saved in jpg format to ", getwd(), "/", folderName, "\n", sep="")
      }
       
      tkinsert(statusWin, "end", statusText)
	  
	  statusText <<- paste("Spatial plot of animal densities saved in jpeg format to ", getwd(), "/", folderName, "\n", sep="")
	  tkinsert(statusWin, "end", statusText)
########### GUI related END ##################
########### Output related START ##############
	  tclarray <- tclArray()
      fname <- paste(folderName, "/summary_stats_", ts,".csv", sep="")
      
      
      if(tclvalue(rbValue78)=="1")  {
        nrowsTclarray <- 8; nrowsTable <- 9 
        write.table(file=fname, row.names=FALSE, col.names=FALSE, sep=",", resTable)
      }else {     #Rashmi
        nrowsTclarray <- 6; nrowsTable <- 7   
        write.table(file=fname, row.names=FALSE, col.names=FALSE, sep=",", resTable[1:3,])  
        write.table(file=fname, append = TRUE, row.names=FALSE, col.names=FALSE, sep=",", resTable[5:7,]) #should be 5:7
        }
      
       # nrowsTclarray <- 6; nrowsTable <- 7   
       # write.table(file=fname, row.names=FALSE, col.names=FALSE, sep=",", resTable[1:7,])  }
        
        
########### Output related END ###############
########### GUI related START ################
      statusText <<- paste("Summary statistics written to ", getwd(), "/", fname, "\n", sep="") 
      tkinsert(statusWin, "end", statusText)
      
      if (Mb == 0)            #Rashmi
      {
      for (i in (0:2))
        for (j in (0:4))
           tclarray[[i,j]] <- resTable[i+1,j+1]
           
      for (i in (4:6))    #should be 4:6
        for (j in (0:4))
           tclarray[[i-1,j]] <- resTable[i+1,j+1]      #i is made to be i-1
        table1 <- tkwidget(two,"table",variable=tclarray,rows=nrowsTable-1,cols=5,titlerows=1,titlecols=1,selectmode="extended",colwidth=20,background="white",state="disabled",borderwidth=2)
      }
      else
      {
      for (i in (0:nrowsTclarray))
        for (j in (0:4))
           tclarray[[i,j]] <- resTable[i+1,j+1]
           table1 <- tkwidget(two,"table",variable=tclarray,rows=nrowsTable,cols=5,titlerows=1,titlecols=1,selectmode="extended",colwidth=20,background="white",state="disabled",borderwidth=2)
      }
         
      

      
            
#      i<-1
#      img <-tkrplot(two, function(){par(bg="white",fin=c(3,3));plot(density(out[,i], na.rm=TRUE),main=resTable[i+1,1],)},hscale=0.75, vscale=0.75)
#      f<-function(...) {
#          val <- as.numeric(tclvalue("i"))
#          if (val!= i) {
#              i <<- val
#             tkrreplot(img)
#          }
#      }
#      s <- tkscale(two, command=f, from=1, to=5, variable="i",showvalue=T, resolution=1, orient="horiz")
      
    		mod1<-""
    		mod2<-""
        if(tclvalue(rbValue34)=="1")
    		{ mod2<-"Spatial Capture-Recapture"
			}else{ mod2<-"Non-Spatial Capture-Recapture" }
    		
        if(tclvalue(rbValue78)=="1")
    		{ mod1<-"Trap response present"
    		}else{ mod1<-"Trap response absent" }
		
		if(tclvalue(rbValue56)=="1")
			{ mod3<-"Half-normal detection function"
			}else{ mod3<-"Negative exponential function" }
    		
    		mod4<-"Bernoulli detection process"
    		
        txt1 <- paste("Input Summary\nArea of each pixel representing a potential home-range center:", tclvalue(acSize), "sq km")
        txt2 <- paste("Model selected: ", mod1, ", ", mod2, ", ", mod3, ", ", mod4, sep="")
        txt3 <- paste("MCMC simulation settings: Iterations -", Global_iterNum, "Burnin -", Global_burn, "Thinning -", Global_skip, "Data Augmentation -", Global_nz)
        lab1 <- tklabel(two, text=txt1, background="white")
        lab2 <- tklabel(two, text=txt2, background="white")
        lab3 <- tklabel(two, text=txt3, background="white")
        tkgrid(lab1, sticky="w")
        tkgrid(lab2, sticky="w")
        tkgrid(lab3, sticky="w")        
                
      tkgrid(table1, sticky="w")
      #tkpack(img,fill="x",pady=4)
      #tkpack(s,fill="x",pady=4)
      tkpack(two,fill="both",expand="yes")

    }  ## End of SCRd.fn 
  
	
    run <-function()
  	{
      if( tclvalue(tkcget(ok1,"-state"))=="normal" | tclvalue(tkcget(ok2,"-state"))=="normal" | tclvalue(tkcget(ok3,"-state"))=="normal") {
          tkmessageBox(message="Please complete input data, model definition and MCMC simulation settings before starting the analysis!",icon="error",type="ok")
          statusText <<- "Please complete input data, model definition and MCMC simulation settings before starting the analysis\n"
          tkinsert(statusWin, "end", statusText) 
          tkfocus(tt)
        }else{
          tkdestroy(two)
          two <<- tkframe(tt, width = 800, height = 300, bg= "white",relief="groove",borderwidth=2)
         	tkpack(two, side = "top")
          
      		ni <-as.integer(Global_iterNum)
      		burn <- as.integer(Global_burn)
      		skip <- as.integer(Global_skip)
      		nz <- as.integer(Global_nz)
      		
      		mod1<-""
      		mod2<-""
			mod3<-""
          if(tclvalue(rbValue34)=="1")
      		{ bsigma=1; mod2<-"Spatial Capture-Recapture"
          }else{ bsigma=0; mod2<-"Non-Spatial Capture-Recapture" }
      		
          if(tclvalue(rbValue78)=="1")
      		{Mb=1; mod1<-"Trap response present"
      		}else{ Mb=0; mod1<-"Trap response absent" }
		
		  if(tclvalue(rbValue56)=="1")
      		{dexp=2; mod3<-"Half-normal detection function"
      		}else{ dexp=1; mod3<-"Negative exponential detection function" }
      		
      		mod4<-"Bernoulli detection process"
      		
          statusText <<- paste("Input Summary\n------------\nArea of each pixel representing a potential home-range center:", tclvalue(acSize), "sq km")
          statusText <<- paste(statusText, "\nModel selected: ", mod1, ", ", mod2, ", ", mod3, ", ", mod4, sep="")
          statusText <<- paste(statusText, "\nMCMC simulation settings: Iterations -", ni, "Burnin -", burn, "Thinning -", skip, "Data Augmentation -", nz)
          tkinsert(statusWin, "end", statusText) 

      		SCRd.fn(ni, burn, skip, bsigma, Mb, nz, dexp)  }
  	}

    stopit <- function()      {
      shouldIStop <<- TRUE    }
 
    exit <- function()       {
      tkdestroy(tt)          }
    
    helpme <- function()
    {
      # Specify the location of the help file and invoke it using the 'exec' command.
      # The assumption is that help file - be it in PDF, HTML, PS etc format - has its
      # "HANDLER" registered in the OS.
      # In simple english, it means on the computer where this R program is supposed to run,
      # an application which can open the help file, has been installed and the OS now knows
      # to use that application to open this help file.

      r_helpfile_location <- paste(.Library, "/SPACECAP/doc/SPACECAP_Manual.pdf", sep="")
      if(!is.null(try(shell.exec(r_helpfile_location), TRUE)))
        tkmessageBox(message=geterrmessage(),icon="error",type="ok") 
    }
   
    tkwm.title(tt,"SPACECAP Ver 1.0.5")
  
  	topMenu <- tkmenu(tt)
  	tkconfigure(tt, menu=topMenu)
  	fileMenu <- tkmenu(topMenu, tearoff=FALSE)
  
  	tkadd(topMenu, "command", label="Run", command=run)
  	tkadd(topMenu, "command", label="Stop", command=stopit)
  	tkadd(topMenu, "command", label="Exit", command=exit)
  	tkadd(topMenu, "command", label="Help", command=helpme)
  	      
    # fileFrame contents
    head1 <- tklabel(fileFrame, text="Input Data", background="light blue", pady=5)
  	tkgrid(head1)
  	  
   	lab1 <- tklabel(fileFrame, text="Select animal capture details file", background="light blue")
   	capt <-tkentry(fileFrame, background="white", width=40, state="disable")
   	button1 <- tkbutton(fileFrame, text="Browse", command=readFile1, width=10)
      
   	lab2 <- tklabel(fileFrame, text="Select trap deployment details file", background="light blue")
   	trap <-tkentry(fileFrame, background="white", width=40, state="disable")
   	button2 <- tkbutton(fileFrame, text="Browse", command=readFile2, width=10)

    lab3 <- tklabel(fileFrame, text="Select potential home-range centers data file", background="light blue")
   	grids <-tkentry(fileFrame, background="white", width=40, state="disable")
   	button3 <- tkbutton(fileFrame, text="Browse", command=readFile3, width=10)

    tkgrid(lab3, sticky="w", padx=5)
  	tkgrid(grids, button3, sticky="w", padx=5)
    tkgrid(lab2, sticky="w", padx=5)
  	tkgrid(trap, button2, sticky="w", padx=5)
    tkgrid(lab1, sticky="w", padx=5)
  	tkgrid(capt, button1, sticky="w", padx=5)

    ac_label1 <- tklabel(fileFrame, text="Specify the area of each pixel (in sq. km)", background="light blue")
    ac_label2 <- tklabel(fileFrame, text="that represents a potential home-range center", background="light blue")
  	acSize <- Global_acSize
    ac_size <-tkentry(fileFrame, textvariable = acSize, background="white", width=10)
  	tkgrid(ac_label1, sticky="w", padx=5)
  	tkgrid(ac_label2, sticky="w", padx=5)
  	tkgrid(ac_size, sticky="w", padx=5)
  	
  	OnOk1 <- function() 
    {
        ## Error checks for input files incorporated here
        ## 1. data is numeric with no NA values
        ## 2. Number of cols and column order is correct
        ## 3. captures cross checked with traps
        errorFlag = 0
        if(is.null(locidso) | is.null(locso) | is.null(grid2500)) {
          print("Error - Input files not selected")
          statusText <<- "Error - Input files not selected\n" 
          tkinsert(statusWin, "end", statusText)                  
          tkmessageBox(message="Error - Input files not selected!",icon="error",type="ok")  
          errorFlag = -1 }
  
        ##Check locidso
        if(dim(locidso)[1]==0)
        {   
          print("Error in animal capture details file - No data or bad file format, csv file expected")
          statusText <<- "Error in animal capture details file - No data or bad file format, csv file expected\n" 
          tkinsert(statusWin, "end", statusText)          
          tkmessageBox(message="Error in animal capture details file - No data or bad file format, csv file expected!",icon="error",type="ok")   
          errorFlag = -1  
        }        
        if(dim(locidso)[2]==3)
        {
          if (is.integer(locidso[,1]) & is.integer(locidso[,2]) & is.integer(locidso[,3]) & sum(is.na(locidso))==0)
          {
            if(toupper(names(locidso)[1])!="LOC_ID" | toupper(names(locidso)[2])!="ANIMAL_ID" | toupper(names(locidso)[3])!="SO") {
               print("Error in animal capture details file - incorrect column sequence or header")
               statusText <<- "Error in animal capture details file - incorrect column sequence or header\n" 
               tkinsert(statusWin, "end", statusText)          
               tkmessageBox(message="Error in animal capture details file - incorrect column sequence or header!",icon="error",type="ok")   
               errorFlag = -1  }
          }else {
            print("Error in animal capture details file - non-integer or missing values")  
            statusText <<- "Error in animal capture details file - non-integer or missing values\n" 
            tkinsert(statusWin, "end", statusText)    
            tkmessageBox(message="Error in animal capture details file - non-integer or missing values!",icon="error",type="ok")  
            errorFlag = -1 }
        }else {
          print("Error in animal capture details file - incorrect number of columns")
          statusText <<- "Error in animal capture details file - incorrect number of columns\n" 
          tkinsert(statusWin, "end", statusText)   
          tkmessageBox(message="Error in animal capture details file - incorrect number of columns!",icon="error",type="ok")  
          errorFlag = -1 }
        ### Check locso
        if(dim(locso)[1]==0)
        {   
          print("Error in trap deployment details file - No data or bad file format, csv file expected")
          statusText <<- "Error in trap deployment details file - No data or bad file format, csv file expected\n" 
          tkinsert(statusWin, "end", statusText)          
          tkmessageBox(message="Error in trap deployment details file - No data or bad file format, csv file expected!",icon="error",type="ok")   
          errorFlag = -1  
        }
        if(is.integer(locso[,1]) & is.numeric(locso[,2]) & is.numeric(locso[,3]) & sum(is.na(locso))==0)
        {
          if(toupper(names(locso)[1])!="LOC_ID" | toupper(names(locso)[2])!="X_COORD" | toupper(names(locso)[3])!="Y_COORD") {
            print("Error in trap deployment details file - incorrect column sequence or header")
            statusText <<- "Error in trap deployment details file - incorrect column sequence or header\n" 
            tkinsert(statusWin, "end", statusText)   
            tkmessageBox(message="Error in trap deployment details file - incorrect column sequence or header!",icon="error",type="ok")  
            errorFlag = -1 }
        }else {
          print("Error in trap deployment details file - non-numeric or missing values")
          statusText <<- "Error in trap deployment details file - non-numeric or missing values\n" 
          tkinsert(statusWin, "end", statusText)   
          tkmessageBox(message="Error in trap deployment details file - non-numeric or missing values!",icon="error",type="ok") 
          errorFlag = -1 }
        nso = (dim(locso)[2]) - 3
        nloc = dim(locso)[1]
        for(i in 4:(nso+3)) { 
          if (sum(locso[,i]>=0)<nloc | sum(locso[,i]<=1)<nloc) {
            statusText <<- paste("Error in trap deployment details file - incorrect trap status in column ", i,", should be 0 or 1\n", sep="") 
            print(statusText)
            tkinsert(statusWin, "end", statusText)   
            tkmessageBox(message=statusText,icon="error",type="ok") 
            errorFlag = -1 } }
        ### Check grids
        if(dim(grid2500)[1]==0)
        {   
          print("Error in potential home range center data file - No data or bad file format, csv file expected")
          statusText <<- "Error in potential home range center data file - No data or bad file format, csv file expected\n" 
          tkinsert(statusWin, "end", statusText)          
          tkmessageBox(message="Error in potential home range center data file - No data or bad file format, csv file expected!",icon="error",type="ok")   
          errorFlag = -1  
        }
        if(dim(grid2500)[2]==3)
        {
          if(is.numeric(grid2500[,1]) & is.numeric(grid2500[,2]) & is.integer(grid2500[,3]) & sum(is.na(grid2500))==0)
          {
            if(toupper(names(grid2500)[1])!="X_COORD" | toupper(names(grid2500)[2])!="Y_COORD" | toupper(names(grid2500)[3])!="HABITAT")  {
               print("Error in potential home range center data file - incorrect column sequence or header")
               statusText <<- "Error in potential home range center data file - incorrect column sequence or header\n" 
               tkinsert(statusWin, "end", statusText)   
               tkmessageBox(message="Error in potential home range center data file - incorrect column sequence or header!",icon="error",type="ok")  
               errorFlag = -1  }
          }else {
            print("Error in potential home range center data file - non-numeric or missing values")
            statusText <<- "Error in potential home range center data file - non-numeric or missing values\n" 
            tkinsert(statusWin, "end", statusText)   
            tkmessageBox(message="Error in potential home range center data file - non-numeric or missing values!",icon="error",type="ok")  
            errorFlag = -1   }
        }else  {
          print("Error in potential home range center data file - incorrect number of columns")
          statusText <<- "Error in potential home range center data file - incorrect number of columns\n" 
          tkinsert(statusWin, "end", statusText)   
          tkmessageBox(message="Error in potential home range center data file - incorrect number of columns!",icon="error",type="ok")  
          errorFlag = -1   }

        #### Check captures vs traps
        len = length(locidso[,1])
        for (i in 1:len)
        {
          loc = locidso[i,1]
          so = locidso[i,3]
          
          if(locso[loc,so+3]==0) {
            cat("Error - mismatch in animal capture details and trap deployment details files : location id", loc, "not deployed on SO", so, fill=TRUE)
            statusText <<- paste("Error - mismatch in animal capture details and trap deployment details files : location id", loc, "not deployed on SO", so)
            tkinsert(statusWin, "end", statusText)   
            tkmessageBox(message=paste("Error - mismatch in animal capture details and trap deployment details files : location id", loc, "not deployed on SO", so),icon="error",type="ok")  
            errorFlag = -1   }
        }    

        Global_acSize <<- as.numeric(tclvalue(acSize))
      	  
        if (is.na(Global_acSize)) {                                              
      	  tkmessageBox(message="Error - potential home-range center area should be numeric!",icon="error",type="ok")
      		tkfocus(ac_size)    
          statusText <<- "Error - potential home-range center area should be numeric\n" 
          tkinsert(statusWin, "end", statusText)  
          errorFlag = -1   } 
      
      ##########
      if(errorFlag==-1) {
        statusText <<- "Input data could not be read\n"
        tkinsert(statusWin, "end", statusText)
      }else {
        statusText <<- "Input data read successfully\n"
        tkinsert(statusWin, "end", statusText)  
        tkconfigure(button1,state="disable")
        tkconfigure(button2,state="disable")
        tkconfigure(button3,state="disable")  
        tkconfigure(ac_size,state="disable")  
        tkconfigure(ok1,state="disable")  }

        #tkentryconfigure(topMenu,1,state="active")
  	}
    
    OnReset1 <- function() 
    {
      tkconfigure(button1,state="active")
      tkconfigure(button2,state="active")
      tkconfigure(button3,state="active")
      tkconfigure(ac_size,state="normal")  
      tkconfigure(ok1,state="active")
    }
    
  	ok1 <-tkbutton(fileFrame,text="   OK   ",command = OnOk1, width=8)
  	reset1 <-tkbutton(fileFrame,text="  Edit  ",command = OnReset1, width=8)
  	tkgrid(ok1, reset1, padx=5,pady=40)
  	tkgrid.configure(ok1, sticky="e")
  	tkgrid.configure(reset1, sticky="w")
  	
    #choiceFrame contents
    head2 <- tklabel(choiceFrame, text="Model Definition", background="light blue", pady=5)
  	tkgrid(head2)
  	  
    rb78_label <- tklabel(choiceFrame,text="Trap response", background = "light blue", padx=10, pady=1)
  	rb7 <- tkradiobutton(choiceFrame, text = "Trap response present", background = "light blue", padx=10, pady=1)
  	rb8 <- tkradiobutton(choiceFrame, text = "Trap response absent", background = "light blue", padx=10, pady=1)
   	rbValue78 <- tclVar("0")
   	tkconfigure(rb7,variable=rbValue78,value="1")
  	tkconfigure(rb8,variable=rbValue78,value="0")
    tkgrid( rb78_label,sticky="w" )
  	tkgrid( rb7,sticky="w" )
  	tkgrid( rb8,sticky="w" )
  
  	rb34_label <- tklabel(choiceFrame,text="Capture-Recapture model", background = "light blue", padx=10, pady=1)
  	rb3 <- tkradiobutton(choiceFrame, text = "Spatial Capture-Recapture", background = "light blue", padx=10, pady=1)
  	rb4 <- tkradiobutton(choiceFrame, text = "Non-Spatial Capture-Recapture", background = "light blue", padx=10, pady=1)
   	rbValue34 <- tclVar("1")
  	tkconfigure(rb3,variable=rbValue34,value="1")
  	tkconfigure(rb4,variable=rbValue34,value="0")
    tkgrid( rb34_label,sticky="w" )
  	tkgrid( rb3,sticky="w" )
  	tkgrid( rb4,sticky="w" )
  
    rb56_label <- tklabel(choiceFrame,text="Detection function", background = "light blue", padx=10, pady=1)
  	rb5 <- tkradiobutton(choiceFrame, text = "Half Normal", background = "light blue",padx=10, pady=1)
  	rb6 <- tkradiobutton(choiceFrame, text = "Negative Exponential", background = "light blue",padx=10, pady=1)
  	rbValue56 <- tclVar("1")
  	tkconfigure(rb5,variable=rbValue56,value="1")
  	tkconfigure(rb6,variable=rbValue56,value="0")
  	tkgrid( rb56_label,sticky="w" )
  	tkgrid( rb5,sticky="w")
  	tkgrid( rb6,sticky="w")
  
    rb12_label <- tklabel(choiceFrame,text="Capture encounters", background = "light blue", padx=10, pady=1)
  	rb1 <- tkradiobutton(choiceFrame, text = "Bernoulli process", background = "light blue", padx=10, pady=1)
  	rb2 <- tkradiobutton(choiceFrame, text = "Poisson process", background = "light blue", padx=10, pady=1)
    rbValue12 <- tclVar("1")
    tkconfigure(rb1,variable=rbValue12,value="1",state="disabled")
  	tkconfigure(rb2,variable=rbValue12,value="0",state="disabled")
  	tkgrid( rb12_label,sticky="w" )
  	tkgrid( rb1,sticky="w" )
  	tkgrid( rb2,sticky="w")

    onOk2 <- function()
    {
      tkconfigure(rb3,state="disable")
      tkconfigure(rb4,state="disable")
	  tkconfigure(rb5,state="disable")
	  tkconfigure(rb6,state="disable")
      tkconfigure(rb7,state="disable")
      tkconfigure(rb8,state="disable")
      tkconfigure(ok2,state="disable")
      statusText <<- "Model definition complete\n"
      tkinsert(statusWin, "end", statusText)
    }
    onReset2 <- function()
    {
      tkconfigure(rb3,state="normal")
      tkconfigure(rb4,state="normal")
	  tkconfigure(rb5,state="normal")
	  tkconfigure(rb6,state="normal")
      tkconfigure(rb7,state="normal")
      tkconfigure(rb8,state="normal")
      tkconfigure(ok2,state="active")
    }
    ok2 <-tkbutton(choiceFrame,text="     OK      ",command=onOk2, width=8)
  	reset2 <-tkbutton(choiceFrame,text="  Edit  ",command=onReset2, width=8)
  	tkgrid(ok2, reset2, padx=5,pady=10)
  	#tkgrid(reset2, padx=10,pady=10, sticky="w")
   	#tkgrid(ok1, reset1, padx=5,pady=40)
  	tkgrid.configure(ok2, sticky="e")
  	tkgrid.configure(reset2, sticky="w")
 
  	
  	#valueFrame
  	head3 <- tklabel(valueFrame, text="MCMC simulations settings", background="light blue", pady=5)
  	tkgrid(head3)
  	
  	rb_7_label1 <- tklabel(valueFrame, text="Specify number of MCMC iterations", background="light blue", padx=5, pady=1)
  	rb_7_label2 <- tklabel(valueFrame, text="(usually a large number [>50000])", background="light blue", padx=5, pady=1)
  	ItrNum <- Global_iterNum
  	Iterations <-tkentry(valueFrame, textvariable = ItrNum, background="white", width=10)
  	tkgrid( rb_7_label1,sticky="w")
    tkgrid( rb_7_label2,sticky="w")  	
  	tkgrid(Iterations,sticky="w", padx=5)
  
  	rb_8_label1 <- tklabel(valueFrame, text="Specify the burn-in period", background="light blue", padx=5, pady=1)
  	rb_8_label2 <- tklabel(valueFrame, text="(no. of initial iterations to be discarded)", background="light blue", padx=5, pady=1)
  	burnIn <- Global_burn
  	BurnIn <-tkentry(valueFrame, textvariable = burnIn, background="white", width=10)
  	tkgrid(rb_8_label1,sticky="w")
  	tkgrid(rb_8_label2,sticky="w")
  	tkgrid( BurnIn,sticky="w", padx=5)
  
  
  	rb_9_label1 <- tklabel(valueFrame, text="Specify thinning rate", background="light blue", padx=5, pady=1)
  	rb_9_label2 <- tklabel(valueFrame, text="(no. of iterations to be skipped when", background="light blue", padx=5, pady=1)
  	rb_9_label3 <- tklabel(valueFrame, text="reporting summary statistics)", background="light blue", padx=5, pady=1)
  	skip <- Global_skip
  	Skip <-tkentry(valueFrame, textvariable = skip, background="white", width=10)
  	tkgrid( rb_9_label1,sticky="w")
  	tkgrid( rb_9_label2,sticky="w")
  	tkgrid( rb_9_label3,sticky="w")
  	tkgrid( Skip,sticky="w" , padx=5)
  
  
  	rb_10_label1 <- tklabel(valueFrame, text="Specify data augmentation value", background="light blue", padx=5, pady=1)
  	rb_10_label2 <- tklabel(valueFrame, text="(an upper limit to the number of individuals", background="light blue", padx=5, pady=1)
  	rb_10_label3 <- tklabel(valueFrame, text="in the area)", background="light blue", padx=5, pady=1)
  	nz <- Global_nz
  	NZ <-tkentry(valueFrame, textvariable = nz, background="white", width=10)
  	tkgrid( rb_10_label1,sticky="w")
  	tkgrid( rb_10_label2,sticky="w")
  	tkgrid( rb_10_label3,sticky="w")
  	tkgrid(NZ,sticky="w", padx=5)
  	#
    tkpack(one,two,three,expand=TRUE, fill="both")
    tkpack(fileFrame,choiceFrame,valueFrame,expand=TRUE, fill="both")

  	OnOk3 <- function() 
    {
      errorFlag <- 0
      Global_iterNum <<- as.numeric(tclvalue(ItrNum))
      if (is.na(Global_iterNum)) {
      	  tkmessageBox(message="Number of iterations should be an integer!",icon="error",type="ok")
      		tkfocus(Iterations)    
          statusText <<- "Error - Number of iterations should be numeric\n" 
          tkinsert(statusWin, "end", statusText)  
          errorFlag <- -1        }
  
      Global_burn <<- as.numeric(tclvalue(burnIn))
      if (is.na(Global_burn))   {
    		 tkmessageBox(message="Burn-in value should be an integer!",icon="error",type="ok")
    		 tkfocus(BurnIn)        
         statusText <<- "Error - Burn-in value should be an integer\n" 
         tkinsert(statusWin, "end", statusText)  
         errorFlag <- -1       }
  
      Global_skip <<- as.numeric(tclvalue(skip))
      if (is.na(Global_skip))   {
         tkmessageBox(message="Thinning rate should be an integer!",icon="error",type="ok")
         tkfocus(Skip)      	  
         statusText <<- "Error - Thinning rate should be an integer\n" 
         tkinsert(statusWin, "end", statusText)  
         errorFlag <- -1        }
  
      Global_nz <<- as.numeric(tclvalue(nz))
      if (is.na(Global_nz))  	  {
  		   tkmessageBox(message="Data augmentation value should be an integer!",icon="error",type="ok")
  		   tkfocus(NZ)  	        
         statusText <<- "Error - Data augmentation value should be an integer\n" 
         tkinsert(statusWin, "end", statusText)  
         errorFlag <- -1        }
         
      if(Global_burn >= Global_iterNum) {
         tkmessageBox(message="Burn-in value should be less than number of iterations!",icon="error",type="ok")      
         tkfocus(BurnIn)        
         statusText <<- "Error - Burn-in value should be less than number of iterations\n" 
         tkinsert(statusWin, "end", statusText)
         errorFlag <- -1               }
         
      if(Global_skip >= Global_iterNum) {
         tkmessageBox(message="Thinning should be less than number of iterations!",icon="error",type="ok")      
         tkfocus(BurnIn)        
         statusText <<- "Error - Thinning should be less than number of iterations\n" 
         tkinsert(statusWin, "end", statusText)
         errorFlag <- -1               }
  
      if(errorFlag==-1) {
        statusText <<- "MCMC simulation settings incorrect\n"
        tkinsert(statusWin, "end", statusText)
      }else {
        tkconfigure(Iterations, state="disabled")
        tkconfigure(BurnIn, state="disabled")
        tkconfigure(Skip, state="disabled")
        tkconfigure(NZ, state="disabled")
        tkconfigure(ok3,state="disabled")   
        statusText <<- "MCMC simulation settings complete\n"
        tkinsert(statusWin, "end", statusText)  }
    }

    onReset3 <- function()
    {
      tkconfigure(Iterations, state="normal")
      tkconfigure(BurnIn, state="normal")
      tkconfigure(Skip, state="normal")
      tkconfigure(NZ, state="normal")
      tkconfigure(ok3,state="active")
      #tkentryconfigure(topMenu,1,state="disabled")
    }
    ok3 <- tkbutton(valueFrame,text="     OK      ",command=OnOk3, width=8)
  	reset3 <-tkbutton(valueFrame,text="  Edit  ",command=onReset3, width=8)
  	tkgrid(ok3, reset3, padx=10,pady=10)
  	#tkgrid(reset3, padx=10,pady=10, sticky="w")
   	tkgrid.configure(ok3, sticky="e")
  	tkgrid.configure(reset3, sticky="w")
############################ GUI related END ################################
   }
