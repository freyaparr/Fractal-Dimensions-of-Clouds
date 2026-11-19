# Fractal Dimenisions of Clouds

If Julia isn't actuall installed 
(Check by typing "julia --version" in terminal):
In Terminal
Run: curl -fsSL https://install.julialang.org | sh
Then run: source ~/.bashrc
And check again: julia --version
RUN julia 

For the Packages required in the code:
Run: julia -e 'using Pkg; Pkg.add(["NCDatasets", "FractalDimensions", "ClimateBase", "ComplexityMeasures", "CairoMakie", "Distances"])'

Finally to run a file:
First run: cd "FOLDERNAME THE FILE IS IN"
e.g cd "Radiative Scaling"
Then run: julia "filename"
e.g julia "Radiative Scaling.jl"

MAKE SURE TO COMMIT AFTER SESSION
To commit:
git add .
git commit -m "description of changes"
git push origin main
EASIER METHOD: 
You can also go to the little graph icon on the left and add the files you want to commit, type in the commit note and then click the button saying commit. (If you have uncommited changes it should come up with a little number in the bottom right corner of the icon). Sync changes after.