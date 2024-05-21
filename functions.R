library(raster)
library(dequer)
library(agrmt)

cloudFree <- function(pixelQARaster)
{
  LS9_clearList <- c(21888,21952,21824)
  clearPixels <- pixelQARaster%in%LS9_clearList
  n_clearPixels <- sum(values(clearPixels),na.rm=TRUE)
  n_allPixels <- sum(!is.na(values(pixelQARaster)))  
  return (n_clearPixels/n_allPixels)
}

getArea <- function(img,shp,imp,impMinLevel=0,withImp=TRUE)
{
  
  # cut the satellite to area
  cropped_img <- crop(img, shp)
  masked_img <- mask(cropped_img,shp)  
  
  if (withImp)
  {
    imp_mask <- imp
    imp_mask[imp_mask>impMinLevel & imp_mask<=100] <- 1
    imp_mask[imp_mask!=1] <- NA
    
    origin(imp_mask) <- origin(masked_img)
    
    masked_img <- masked_img * imp_mask
  }
  
  return (masked_img)
}

extractArea <- function(wd,areaMask,impfile,cloudThreshold=0.95)
{
  setwd(wd)
  files <- list.files(path="satImages/",pattern="^LC0.*B10.TIF$")
  nr <- length(files)
  
  shape_mask <- shapefile(paste("limit/",areaMask,sep=""))
  
  imp <- raster(impfile) 
  cropped_imp <- crop(imp, shape_mask)
  masked_imp <- mask(cropped_imp,shape_mask)

  for(i in 1:nr)
  {
    print(paste(i,"/",nr,sep=""))
    
    baseName <- substr(files[i],1,nchar(files[i])-10)
    QAName <- paste(baseName,"QA_PIXEL.TIF",sep="")
    QA <- raster(paste("satImages/",QAName,sep=""))
    
    cropped_QA <- crop(QA, shape_mask)
    masked_QA <- mask(cropped_QA,shape_mask)
    
    cl <- cloudFree(masked_QA)
    cl <- round(cl,3)
    
    if (cl>cloudThreshold)
    {
      
      B10 <- raster(paste("satImages/",files[i],sep=""))
      actualImg <- getArea(B10,shape_mask,masked_imp)
      LST <- actualImg * 0.00341802 + 149 - 273.15
      writeRaster(LST,paste("LST/LST_IMPYES_",files[i],sep=""),overwrite=TRUE)
      
      actualImg <- getArea(B10,shape_mask,masked_imp,withImp = FALSE)
      LST <- actualImg * 0.00341802 + 149 - 273.15
      writeRaster(LST,paste("LST/LST_IMPNO_",files[i],sep=""),overwrite=TRUE)
    }
  }
}

getHotSpots <- function(wd,fn,no=10,minAcceptedValueTh=95,minThreshold=98)
{
  r <- raster(paste(wd,"/LST/",fn,".tif",sep="")) 
  rm <- as.matrix(r)
  visited <- rm 
  visited[,] <- 0  
  dimR <- nrow(rm) 
  dimC <- ncol(rm) 
  v=getValues(r)
  qMin <- quantile(v,minAcceptedValueTh/100,na.rm=TRUE) 
  minThValue <- quantile(v,minThreshold/100,na.rm=TRUE) 
  for (p in 1:no) {    
    v <- getValues(r)
    rm <- as.matrix(r) 
    ma <- max(v,na.rm=TRUE)
    if (ma > qMin) {  
      pma <- which(v==ma)
      row <- pma[1]%/%dimC + 1 
      col <- pma[1]%%dimC
      if (col == 0) {
        row <- row - 1
        col <- dimC
      }
      qu <- stack() 
      push(qu,c(row,col))
      visited[row,col] <- p       
      while (length(qu)>0) {
        head <- pop(qu) 
        ro <- head[1] 
        co <- head[2] 
        for (i in -1:1)      
          for (j in -1:1)
            if (i!=j & (i == 0 | j==0 ))
              if ((ro+i)>=1 & (ro+i)<=dimR & (co+j)>=1 & (co+j)<=dimC)
                if (!is.na(visited[ro+i,co+j]) & !is.na(rm[ro+i,co+j]))
                  if ( visited[ro+i,co+j]==0 & rm[ro+i,co+j]>qMin) {
                    th <- sum((visited == p)*rm,na.rm=TRUE) + rm[ro+i,co+j]
                    th <- th / (sum(visited == p) + 1) 
                    if (th > minThValue) {
                      push(qu,c(ro+i,co+j)) 
                      visited[ro+i,co+j] <- p 
                    }
                  }
      }
      r <- setValues(r, rm * (visited == 0))
    }
  }
  rnew <- setValues(r,visited)
  writeRaster(rnew,paste(wd,"/hotSpots/HS_",no,"_",minAcceptedValueTh,"_",minThreshold,"-",fn,sep=""),overwrite=TRUE)
  return (rnew)
}

getCoolSpots <- function(wd,fn,no=10,minAcceptedValueTh=95,minThreshold=98)
{
  r <- raster(paste(wd,"/LST/",fn,".tif",sep="")) 
  rm <- as.matrix(r)
  visited <- rm 
  visited[,] <- 0 
  dimR <- nrow(rm) 
  dimC <- ncol(rm) 
  v=getValues(r)
  qMin <- quantile(v,1-minAcceptedValueTh/100,na.rm=TRUE)
  minThValue <- quantile(v,1-minThreshold/100,na.rm=TRUE)
  for (p in 1:no)
  {
    v <- getValues(r)
    rm <- as.matrix(r)
    mi <- minnz(na.omit(v))
    if (mi < qMin & mi < minThValue)
    {  
      pma <- which(v==mi)
      row <- pma[1]%/%dimC + 1
      col <- pma[1]%%dimC
      if (col == 0)
      {
        row <- row - 1
        col <- dimC
      }
      qu <- stack()
      push(qu,c(row,col))
      visited[row,col] <- p
      while (length(qu)>0)
      {
        head <- pop(qu)
        ro <- head[1]
        co <- head[2]
        for (i in -1:1)
          for (j in -1:1)
            if (i!=j & (i == 0 | j==0 ))
              if ((ro+i)>=1 & (ro+i)<=dimR & (co+j)>=1 & (co+j)<=dimC)
                if (!is.na(visited[ro+i,co+j]) & !is.na(rm[ro+i,co+j]))
                  if ( visited[ro+i,co+j]==0 & rm[ro+i,co+j]<qMin)
                  {
                    th <- sum((visited == p)*rm,na.rm=TRUE) + rm[ro+i,co+j]
                    th <- th / (sum(visited == p) + 1)
                    if (th < minThValue)
                    {
                      push(qu,c(ro+i,co+j))
                      visited[ro+i,co+j] <- p
                    }
                  }
      }
 
      r <- setValues(r, rm * (visited == 0))
    }
  }
  rnew <- setValues(r,visited)
  writeRaster(rnew,paste(wd,"/coolSpots/CS_",no,"_",minAcceptedValueTh,"_",minThreshold,"-",fn,sep=""),overwrite=TRUE)
  return (rnew)
}

hotspotPersistence <- function(wd,type,limit=10,partNo=3)
{
  setwd(wd) 
  if (type=="HS") 
    files <- c(list.files(path="hotSpots/",pattern=paste("HS_",limit,".*.IMPNO.*.grd$",sep=""))) 
  else 
    if (type=="IMP_HS") 
      files <- c(list.files(path="hotSpots/",pattern=paste("HS_",limit,".*.IMPYES.*.grd$",sep="")))
    
    nr <- length(files)
    r <- raster(paste("hotSpots/",files[1],sep=""))  
    rez <- r 
    
    rez[rez>limit] <- 0 
    rez[rez>0] <- 1 
    
    for(i in 2: nr) 
    { 
      r <- raster(paste("hotSpots/",files[i],sep=""))  
      r[r>limit] <- 0 
      r[r>0] <- 1 
      rez <- rez + r 
    } 
    rez <- rez / nr 
    
    rez[rez==0] <- NA 
    
    limits <- c(0) 
    stp <- 1/partNo 
    for(i in 1:(partNo-1)) 
      limits <- c(limits,i*stp) 
    limits <- c(limits,1) 
    
    for(i in partNo:1) 
      rez[rez>limits[i] & rez<=limits[i+1]] <- i 
    
    writeRaster(rez,paste("persistence/T_nr",nr,"_limit",limit,"_part",partNo,"_",
                          type,sep=""), overwrite=TRUE)
}

coolspotPersistence <- function(wd,type,limit=10,partNo=3)
{
  setwd(wd) 
  if (type=="CS") 
    files <- c(list.files(path="coolSpots/",pattern=paste("CS_",limit,".*.IMPNO.*.grd$",sep=""))) 
  else 
    if (type=="IMP_CS") 
      files <- c(list.files(path="coolSpots/",pattern=paste("CS_",limit,".*.IMPYES.*.grd$",sep="")))
    
    nr <- length(files)
    r <- raster(paste("coolSpots/",files[1],sep=""))  
    rez <- r 
    
    rez[rez>limit] <- 0 
    rez[rez>0] <- 1 
    
    for(i in 2: nr) 
    { 
      r <- raster(paste("coolSpots/",files[i],sep=""))  
      r[r>limit] <- 0 
      r[r>0] <- 1 
      rez <- rez + r 
    } 
    rez <- rez / nr 
    
    rez[rez==0] <- NA 
    
    limits <- c(0) 
    stp <- 1/partNo 
    for(i in 1:(partNo-1)) 
      limits <- c(limits,i*stp) 
    limits <- c(limits,1) 
    
    for(i in partNo:1) 
      rez[rez>limits[i] & rez<=limits[i+1]] <- i 
    
    writeRaster(rez,paste("persistence/T_nr",nr,"_limit",limit,"_part",partNo,"_",
                          type,sep=""), overwrite=TRUE)
}

getHighestIntensity <- function(n,no=3,bitNo=6)
{
  fq =c()
  for(i in 0:(no-1))
  {
    fq <- c(fq,bitwAnd(bitwShiftR(n,i*bitNo),2^bitNo-1))
  }
  mx <- 0
  pmx <- 0
  for(i in 1:no)
    if (fq[i]>mx)
    {
      mx <- fq[i]
      pmx <- i
    }
  return (pmx)
}

hotspotOverallIntensity <- function(wd,type,limit=10,partNo=3,bitNo=6)
{
  setwd(wd)
  if (type=="HS")
    files <- c(list.files(path="hotSpots/", paste("HS_",limit,".*.IMPNO.*.grd$",sep="")))
  else
    if (type=="IMP_HS")
      files <- c(list.files(path="hotSpots/", paste("HS_",limit,".*.IMPYES.*.grd$",sep="")))
    nr <- length(files)
    r <- raster(paste("hotSpots/",files[1],sep="")) 
    r[is.na(r)]<-0
    rez <- r
    rez[rez>0] <- 0
    stp <- limit/partNo
    limits <- c(0)
    for(i in 1:(partNo-1))
      limits <- c(limits,round(i*stp))
    limits <- c(limits,limit)
    for(i in 1: nr)
    {
      r <- raster(paste("hotSpots/",files[i],sep="")) 
      r[is.na(r)]<-0
      r[r>limit] <- 0
      for(j in 1:partNo)
      {
        r[r>limits[j] & r<=limits[j+1]] <- 2^((j-1)*bitNo)
      }
      rez <- rez + r
    }
    vs <- getValues(rez)
    for(i in 1: length(vs))
      vs[i] <- getHighestIntensity(vs[i],partNo,bitNo)
    rez <- setValues(rez,vs)
    rez[rez==0] <- NA         
    writeRaster(rez,paste("overall_Intensity/I_nr",nr,"_limit",limit,"_part",partNo,
                          "_",type,sep=""),overwrite=TRUE)
}

coolspotOverallIntensity <- function(wd,type,limit=10,partNo=3,bitNo=6)
{
  setwd(wd)
  if (type=="CS")
    files <- c(list.files(path="coolSpots/", paste("CS_",limit,".*.IMPNO.*.grd$",sep="")))
  else
    if (type=="IMP_CS")
      files <- c(list.files(path="coolSpots/", paste("CS_",limit,".*.IMPYES.*.grd$",sep="")))
    nr <- length(files)
    r <- raster(paste("coolSpots/",files[1],sep="")) 
    r[is.na(r)]<-0
    rez <- r
    rez[rez>0] <- 0
    stp <- limit/partNo
    limits <- c(0)
    for(i in 1:(partNo-1))
      limits <- c(limits,round(i*stp))
    limits <- c(limits,limit)
    for(i in 1: nr)
    {
      r <- raster(paste("coolSpots/",files[i],sep="")) 
      r[is.na(r)]<-0
      r[r>limit] <- 0
      for(j in 1:partNo)
      {
        r[r>limits[j] & r<=limits[j+1]] <- 2^((j-1)*bitNo)
      }
      rez <- rez + r
    }
    vs <- getValues(rez)
    for(i in 1: length(vs))
      vs[i] <- getHighestIntensity(vs[i],partNo,bitNo)
    rez <- setValues(rez,vs)
    rez[rez==0] <- NA         
    writeRaster(rez,paste("overall_Intensity/I_nr",nr,"_limit",limit,"_part",partNo,
                          "_",type,sep=""),overwrite=TRUE)
}

getAllSpots <- function(wd,no=10,minAcceptedValueTh=95,minThreshold=98)
{
  setwd(wd)
  files <- list.files(path="LST/",pattern="^LST.*.C0.*B10.TIF$")  
  nr <- length(files)
  
  for(i in 1:nr)
  {
    getHotSpots(wd,substr(files[i],1,nchar(files[i])-4),no,minAcceptedValueTh,minThreshold)
    getCoolSpots(wd,substr(files[i],1,nchar(files[i])-4),no,minAcceptedValueTh,minThreshold)
  }
}


