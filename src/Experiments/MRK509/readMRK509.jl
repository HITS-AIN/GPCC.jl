function readMRK509(year)

   year == 2016 || year == 2017 ? nothing : error("year must be either 2016 ir 2017")

   stryear = year == 2016 ? "2016" : "2017"

   files = ["Experiments/MRK509/"*stryear*"/Mrk509.4300.6p0.comb.mjy.dat",
            "Experiments/MRK509/"*stryear*"/Mrk509.5700.6p0.comb.mjy.dat",
            "Experiments/MRK509/"*stryear*"/Mrk509.6200.6p0.comb.mjy.dat",
            "Experiments/MRK509/"*stryear*"/Mrk509.7000.6p0.comb.mjy.dat"]

   lambda = [4300.0;
             5700.0;
             6200.0;
             7000.0]


   alldata = map(readdlm, files)


   m = mean(alldata[1][:,1])

   tobs = [(a[:,1].-m) for a in alldata]
   yobs = [a[:,2]      for a in alldata]
   σobs = [a[:,3]      for a in alldata]


   figure(-1); cla()

   clrs = ["b","g","r","m"]
   for i in 1:length(files)
      plot(tobs[i], yobs[i], "o"*clrs[i])
   end

   return lambda, tobs, yobs, σobs

end
