
## converts pixel data into a cloud or not a cloud 
function makeMask(field, threshold; below=true)
    if below
        return field .< threshold #used if temp
    else
        return field .> threshold #used if light
    end
end

 T = [250 200; 210 250]

 makeMask(T,215)

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

function filterClouds(AreaList,perimetreList,touchingEdge; minPixels=24)
    correctClouds = (AreaList .>= minPixels).& (.!touchingEdge)

    return AreaList[correctClouds],perimetreList[correctClouds]
end

using Statistics

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

    return fractalDimension
end