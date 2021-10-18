function readMRK279()

   files = ["Experiments/MRK279/Mrk279.NB4300.6p0.2017.mjy.dat",
            "Experiments/MRK279/Mrk279.NB5700.6p0.2017.mjy.dat",
            "Experiments/MRK279/Mrk279.NB6200.6p0.2017.mjy.dat",
            "Experiments/MRK279/Mrk279.NB7000.6p0.2017.mjy.dat"]

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
