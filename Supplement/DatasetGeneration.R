##########################################
###
### write a function that simulates one dimensional
### data with several levels. The length of the levels is a poisson distribution + 1
### and the blocks will be distributed uniformly
###
### Input parameters:   n: number of data points
###                     values: vector of values of the non-zero data points
###                     numBlocks: vector with the number of blocks per this value
###                     lenBlocks: vector with a lambda value for the poisson distribution +1
###                                 for the length of the blocks
###                     sigma: value with the sd of the added noise
### Output: a vector of length n with the simulated data
###########################################

simulateOneDim = function(n, values, numBlocks, lenBlocks, sigma)
{
    res = rep(0,n)
    
    numValues = length(values)
    for(i in 1:numValues) ### go through all different possible values
    {
        ### insert the blocks
        for(j in 1:numBlocks[i])
        {
            thisBlockLen = rpois(1,lenBlocks[i])+1
            thisBlockPos = sample(n-thisBlockLen+1, 1)
            res[thisBlockPos:(thisBlockPos+thisBlockLen-1)]=values[i]
        }
    }
    return(res+ rnorm(n, sd=sigma))
}


############################################################
###
### write a function that simulates two-dimensional artificial data
### data with several levels. The length of the levels in both dimensions is a poisson 
### distribution + 1
### and the blocks will be distributed uniformly
###
### Input parameters:   nx,ny: number of data points in each dimension
###                     values: vector of values of the non-zero data points
###                     numBlocks: vector with the number of blocks per this value
###                     lenXBlocks: vector with a lambda value for the poisson distribution +1
###                                 for the length of the blocks in the X-dimension
###                     lenYBlocks: same as lenXBlocks for the Y-dimension
###                     sigma: value with the sd of the added noise
### Output: a vector of length n with the simulated data
############################################################

simulateTwoDim = function(nx, ny, values, numBlocks, lenXBlocks, lenYBlocks, sigma)
{
    res = matrix(rep(0,nx*ny), ncol=nx)
    
    numValues = length(values)
    for(i in 1:numValues) ### go through all different possible values
    {
        ### insert the blocks
        for(j in 1:numBlocks[i])
        {
            thisBlockLenX = rpois(1,lenXBlocks[i])+1
            thisBlockLenY = rpois(1,lenYBlocks[i])+1
            thisBlockPosX = sample(nx-thisBlockLenX+1, 1)
            thisBlockPosY = sample(ny-thisBlockLenY+1, 1)
            rangeX = thisBlockPosX:(thisBlockPosX+thisBlockLenX-1)
            rangeY = thisBlockPosY:(thisBlockPosY+thisBlockLenY-1) 
            res[rangeX, rangeY]=values[i]
        }
    }
    return(res+ rnorm(nx*ny, sd=sigma))
}

