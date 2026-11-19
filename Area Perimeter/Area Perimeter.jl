using NCDatasets,CairoMakie #nc datasets is used to get data - moved from python to julia

botanyFile = "/workspaces/Fractal-Dimensions-of-Clouds/dales_radiation_rlut_data.nc"
ds = NCDataset(botanyFile)

## this section splits clouds up so we end up finding the area and perimetre of individual clouds rather than a conglomerate
function splitClouds(mask)
    zArray = zeros(Int, size(mask)) #creates an array equal to image size, of zeroes
    rows, cols = size(mask)
    tempStack = [] #stack for pixels that need to be checked
    totCloudCount = 0

    for i in 1:rows, j in 1:cols #iterate over every pixel
        
        if mask[i, j] && zArray[i, j] == 0 #checks if we have found an unlabelled cloud
            totCloudCount +=1 #this becomes the new cloud num
            zArray[i, j] = totCloudCount
            push!(tempStack, (i,j)) #ensures neighbours are checked

            while !isempty(tempStack)
                tI,tJ = pop!(tempStack) #continues until stack is empty (all nb are checked)

                for (fI,fJ) in ((-1,0),(1,0),(0,-1),(0,1)) #checks nbs
                    nI,nJ = tI+fI,tJ+fJ


                    if 1<= nI <= rows && 1 <= nJ <= cols && mask[nI,nJ] == 1 && zArray[nI,nJ] == 0 #checks if nb should join cloud if is inside gird, nb is cloud and is unlabelled
                        zArray[nI,nJ] = totCloudCount
                        push!(tempStack, (nI,nJ)) #makes sure nbs get checked
                    end

                end
            end
        end
    end

    return zArray, totCloudCount
end



 ## fill in the clouds so that the inside of the cloud doesnt reflect as multiple
function fillClouds(mask)
    rows, cols = size(mask)
    tempStack = []
    reachedPixels = falses(rows,cols)

    for i in 1:rows, j in 1:cols
        if (i == 1 || i == rows || j==1 || j == cols) && !mask[i,j] #checks if a pixel is border or is sky
            reachedPixels[i,j] = true
            push!(tempStack, (i,j)) #allows this pixels neighbours to be checked

            while !isempty(tempStack)
                tI,tJ = pop!(tempStack)

                for (fI,fJ) in ((-1,0),(1,0),(0,-1),(0,1))
                    nI,nJ = tI+fI,tJ+fJ #checks the pixels neighbours


                    if 1<= nI <= rows && 1 <= nJ <= cols && !mask[nI,nJ] == 1 && !reachedPixels[nI,nJ] #checks if nb row or col are in the grid
                        reachedPixels[nI,nJ] = true
                        push!(tempStack, (nI,nJ))
                    end

                end
            end
        end
    end

    return .!reachedPixels
end




## area perimetre


function areaPerim(zArray,totCloudCount)
    AreaList = zeros(Int, totCloudCount)
    perimetreList = zeros(Int, totCloudCount)
    rows, cols = size(zArray)

    touchingEdge = falses(totCloudCount)

##iterate through the rows and collumns, check if they are cloud or sky
    for i in 1:rows, j in 1:cols
        k = zArray[i, j]
        if k == 0 
            nothing
        else
            AreaList[k] +=1

            for (fI,fJ) in ((-1,0),(1,0),(0,-1),(0,1)) #checks nbs
                    nI,nJ = i+fI,j+fJ


                    if !(1<= nI <= rows && 1 <= nJ <= cols)
                        perimetreList[k] += 1
                        touchingEdge[k] = true
                    elseif zArray[nI,nJ] != k
                        perimetreList[k] += 1
                    end
                end
            end
        end
    return AreaList, perimetreList,touchingEdge
end



#cloud filter uses same rules as christen and driver (24 pixels away from edge)

function filterClouds(AreaList, perimetreList, touchingEdge; minPixels = 24)
    correctClouds = (AreaList .>= minPixels) .& (.!touchingEdge)
    return AreaList[correctClouds], perimetreList[correctClouds]
end

function getD(AreaList,perimetreList)
    xList = log10.(AreaList)
    yList = log10.(perimetreList)
    xBar = sum(xList)/length(xList)
    yBar = sum(yList)/length(yList)
    sum1 = 0
    sum2 = 0
    for i in 1:length(xList)
        sum1 += (xList[i]-xBar)*(yList[i]-yBar)
        sum2+= (xList[i]-xBar)*(xList[i]-xBar)
    end
    grad = (sum1)/sum2
    fractalDimension = 2*grad
    intercept = yBar - grad*xBar
    return fractalDimension,intercept
end


fieldArray = Matrix{Matrix{Float64}}(undef, 5, 6)
for i in 1:5
    for j in 1:6
        fieldArray[i,j] = Float64.(ds["lwp"][:, :, j, i])
    end
end

close(ds)


function fracDimCalc(field,threshold)
    cloudMask = makeMask(field,threshold;below=false)
    filledClouds = fillClouds(cloudMask)
    indClouds,n = splitClouds(filledClouds)
    aList,pList,tEdge = areaPerim(indClouds,n)
    aList,pList = filterClouds(aList,pList,tEdge)
    fracDim,intercept = getD(aList,pList)

    return fracDim,intercept,length(aList),aList,pList
end


function plotFracDim(aList,pList,intercept,fracDim)

    aLogs = log10.(aList)
    pLogs = log10.(pList)

    xMin = minimum(aLogs)
    xMax = maximum(aLogs)
    yArr = []
    xArr= []
    interval = (xMax-xMin)/50
    for i in 1:50
        y = intercept + (0.5*fracDim)*xMin
        push!(yArr,y)
        push!(xArr,xMin)
        xMin+=interval
    end

    fig = Figure()
    ax = Axis(fig[1,1],title = "Log Area against Log Perimeter",xlabel = "Log Area",ylabel = "Log Perimeter")
    scatter!(ax,aLogs,pLogs, color = (:green, 0.2), markersize = 12)
    lines!(ax,xArr,yArr, color = (:black, 0.8))
    fig

end

fracDimsFinal = []
interceptFinal = []
totalAreaList = []
totalPerimList = []

for i in 1:5
    for j in 1:6
        D, c, n, aL, pL = fracDimCalc(fieldArray[i,j], 0.01)
        push!(fracDimsFinal,D)
        push!(interceptFinal,c)
        append!(totalAreaList, aL)
        append!(totalPerimList, pL)
    end
end

validFracDims = filter(!isnan, fracDimsFinal)
validInter = filter(!isnan, interceptFinal)

averageDim = sum(validFracDims)/length(validFracDims)
averageInter = sum(validInter)/length(validInter)

poolDim, poolInter = getD(totalAreaList, totalPerimList)
plotFracDim(totalAreaList,totalPerimList,poolInter, poolDim)