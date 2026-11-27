






using NCDatasets,CairoMakie #nc datasets is used to get data - moved from python to julia



## converts pixel data into a cloud or not a cloud - every pixel becomes either sky or cloud
function makeMask(field, threshold; below=true) #takes a 1536*1536 matrix of numbers, threshold refers to which numbers stay
    if below
        return field .< threshold #used if temp
    else
        return field .> threshold #used if light
    end
end


## this section splits clouds up so we end up finding the area and perimetre of individual clouds rather than a conglomerate
function splitClouds(mask)#mask is the filled in clouds from fill clouds
    zArray = zeros(Int, size(mask)) #creates an array equal to image size, of zeroes, used to see if a cloud is reached and to hold which cloud it is
    rows, cols = size(mask)# two variable that describes the size of the image
    tempStack = [] #stack for pixels that need to be checked
    totCloudCount = 0 #the total amount of clouds per scene

    for i in 1:rows, j in 1:cols # a nested loop to iterate over every pixel
        
        if mask[i, j] == true && zArray[i, j] == 0 #checks if [i j]  is both cloud and unvisited
            totCloudCount +=1 #this becomes the new cloud num
            zArray[i, j] = totCloudCount#sets this as cloud 
            push!(tempStack, (i,j)) #ensures neighbours are checked

            while !isempty(tempStack)#iterates over the entire stack list
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

#This code can tell between 2 clouds because it only loops when every neighbour is checked and it has found an unlabeled pixel that is cloud
end



 ## fill in the clouds so that the inside of the cloud doesnt incorrectly add to the clouds perimeter
function fillClouds(mask)
    rows, cols = size(mask) # creates two numbers that describe the size of the image
    tempStack = [] #an emoty array
    reachedPixels = falses(rows,cols) #dimension is rows * cols, where every element = false

    for i in 1:rows, j in 1:cols# creates a nested loop

        if (i == 1 || i == rows || j==1 || j == cols) && !mask[i,j] #checks if a pixel is on the border and if it is sky
            reachedPixels[i,j] = true #turns any pixel that has been checked from false to true
            push!(tempStack, (i,j)) #adds all the reached pixels into a stack

            while !isempty(tempStack)#iterates across every reached pixel
                tI,tJ = pop!(tempStack) #removes it from the stack and stores it

                for (fI,fJ) in ((-1,0),(1,0),(0,-1),(0,1))
                    nI,nJ = tI+fI,tJ+fJ #checks the pixels neighbours, for any pixel above below left or right using 4 o


                    if 1<= nI <= rows && 1 <= nJ <= cols && !mask[nI,nJ] == 1 && !reachedPixels[nI,nJ] #checks if the popped off pixels nb is on the image, is sky, and hasnt been flooded
                        reachedPixels[nI,nJ] = true #turns it to reached, and adds it to the stack
                        push!(tempStack, (nI,nJ))
                    end

                end
            end
        end
    end

    return .!reachedPixels # returns all pixels that werent flooded, which means they are either a cloud pixel or inside a set of cloud. This makes it so that any unfilled in bits become cloud.
end




## area perimeter


function areaPerim(zArray,totCloudCount)
    AreaList = zeros(Int, totCloudCount)#two arrays of length totcloudcount with only 0s, used to gewt area and perim of each cloud
    perimetreList = zeros(Int, totCloudCount)
    rows, cols = size(zArray)

    touchingEdge = falses(totCloudCount)#checks if each cloud is touchign the edge of the image

#iterate through the rows and collumns, check if they are cloud or sky
    for i in 1:rows, j in 1:cols
        k = zArray[i, j]
        if k == 0 #checks if its sky
            nothing
        else
            AreaList[k] +=1 #iterates the area list by 1

            for (fI,fJ) in ((-1,0),(1,0),(0,-1),(0,1)) #checks nbs as previously
                    nI,nJ = i+fI,j+fJ


                    if !(1<= nI <= rows && 1 <= nJ <= cols)#this pixel is on the edge of the image so we increase perimeter and count it as on the edge
                        perimetreList[k] += 1
                        touchingEdge[k] = true
                    elseif zArray[nI,nJ] != k #checks if its not part of the current cloud - either sky or another clouds ooutline, so increases perimeter by one.
                        perimetreList[k] += 1
                    end
                end
            end
        end
    return AreaList, perimetreList,touchingEdge
end



#cloud filter uses same rules as christen and driver (24 pixels away from edge)

function filterClouds(AreaList, perimetreList, touchingEdge; minPixels = 24)
    correctClouds = (AreaList .>= minPixels) .& (.!touchingEdge) #checks that every cloud meets the min area size and isnt touching the edge and returns them
    return AreaList[correctClouds], perimetreList[correctClouds]
end

function getD(AreaList,perimetreList)
    xList = log10.(AreaList)#logs all of area and perim
    yList = log10.(perimetreList)
    xBar = sum(xList)/length(xList)#gets avergaes
    yBar = sum(yList)/length(yList)
    sum1 = 0
    sum2 = 0
    for i in 1:length(xList)
        sum1 += (xList[i]-xBar)*(yList[i]-yBar)
        sum2+= (xList[i]-xBar)*(xList[i]-xBar)##calculates least square slope - represents all points better than gradient

    end
    grad = (sum1)/sum2
    fractalDimension = 2*grad
    intercept = yBar - grad*xBar
    return fractalDimension,intercept
end

using NCDatasets#opens data and makes sure the data is actual datan 
botanyFile = "/Users/barneywelford/Downloads/Day and Night Data/rlut_12am_feb3.nc"
ds = NCDataset(botanyFile)
fileData = ds["cloud_mask"].var



goodData=[]
for i in 1:103
    push!(goodData,i)
end



function fracDimCalc(field,threshold)
    mask = field .==1
    filledClouds = fillClouds(mask)
    indClouds,n = splitClouds(filledClouds)
    aList,pList,tEdge = areaPerim(indClouds,n)
    aList,pList = filterClouds(aList,pList,tEdge)
    fracDim,intercept = getD(aList,pList)

    return fracDim,intercept,length(aList),aList,pList
end#actually runs the data pipeline



function plotFracDim(aList,pList,intercept,fracDim)#plots the data, as well as using the given intercept and fracDim to get a line of best fit


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

# fracDim,intercept,lengthAlist,aList,pList = fracDimCalc(fileData[:,:,6],0.01)
# plotFracDim(aList,pList,intercept,fracDim)


totalFracDim = []
totalIntercept = []
totalN = []
for i in goodData
    fracDim,intercept,lengthAlist,aList,pList = fracDimCalc(fileData[:,:,i],0.01)
    push!(totalFracDim,fracDim)
    push!(totalIntercept,intercept)
    push!(totalN,lengthAlist)
    println(fracDim," ",intercept," ")
end
using Statistics

validFracDims = filter(!isnan, totalFracDim)
validInter = filter(!isnan, totalIntercept)






using CSV,DataFrames

sliceData = DataFrame(simulation = goodData,fractalDimension=totalFracDim,intercepts=totalIntercept)
CSV.write("areaPerimData.csv",sliceData)
out = joinpath(homedir(), "Downloads", "rlut_12am_feb3.nc_FD.csv")
CSV.write(out, sliceData)
