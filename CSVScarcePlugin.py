import PyPluMA

def isnumber(mystring):
    digitcount = 0
    for i in range(0, len(mystring)):
        c = mystring[i]
        if ((not c.isdigit()) and (c != '.') and (not (i == 0 and c == '-'))):
            return False
        elif (c.isdigit()):
            digitcount += 1
    if (digitcount != 0):
       return True


class CSVScarcePlugin:
    def input(self, filename):
      self.parameters = dict()
      paramfile = open(filename, 'r')
      for line in paramfile:
         contents = line.split('\t')
         self.parameters[contents[0]] = contents[1].strip()

      self.otu_file = open(PyPluMA.prefix()+"/"+self.parameters["OTU"], 'r')
      sample_file =open(PyPluMA.prefix()+"/"+self.parameters["META"], 'r')
      self.threshold = float(self.parameters["threshold"])
      #if ("keepifone" in self.parameters):
      #    self.keepifone = True
      #else:
      #    self.keepifone = False
      #Idea: Zero out (don't remove) OTUs in less than 50% of a category set
      # Then once everything is done, remove those with zero over all samples

      self.categories = dict()
      sample_file.readline()
      for line in sample_file:
          contents = line.strip().split(',')
          self.categories[contents[0]] = contents[1]
      
      self.lines = []
      firstline = self.otu_file.readline().strip()
      self.taxa = firstline.split(',')
      self.counts = dict()
      for i in range(1, len(self.taxa)):
          taxon = self.taxa[i]
          self.counts[taxon] = dict()
      
      self.totals = dict()
      
      for line in self.otu_file:
          contents = line.strip().split(',')
          self.lines.append(contents)
          sample = contents[0]
          category = self.categories[sample]
          if category not in self.totals:
              self.totals[category] = 0
          self.totals[category] += 1
          for i in range(1, len(contents)):
              taxon = self.taxa[i]
              if (category not in self.counts[taxon]):
                 self.counts[taxon][category] = 0
              if (not isnumber(contents[i]) or float(contents[i]) != 0):
                  self.counts[taxon][category] += 1
                      
           
    def run(self):
       # Two different possibilities.
       # keepifone = False:
       # Keep taxa in the category they broke the threshold and zero out other places
       # keepifone = True:
       # Do not zero out other places
       # Either way, remove all taxa that never broke the threshold in any category
       self.toRemove = []
       for taxon in self.counts:
           brokeThresh = []
           for category in self.counts[taxon]:
             if (float(self.counts[taxon][category])/self.totals[category] > self.threshold):
                brokeThresh.append(category)
           if (len(brokeThresh) == 0):
              self.toRemove.append(taxon)

    def output(self, filename):
      otu_filter_file = open(filename, 'w')
      finaltaxa = []
      lastpos = -1
      for i in range(len(self.taxa)):
          taxon = self.taxa[i]
          if (taxon not in self.toRemove):
              finaltaxa.append(taxon)
              lastpos = i
      for i in range(len(finaltaxa)):
          otu_filter_file.write(finaltaxa[i])
          if (i != len(finaltaxa)-1):
              otu_filter_file.write(',')
          else:
              otu_filter_file.write('\n')
      for line in self.lines:
          for i in range(lastpos+1):
              taxon = self.taxa[i]
              if (taxon not in self.toRemove):
                  otu_filter_file.write(line[i])
                  if (i != lastpos):
                      otu_filter_file.write(',')
                  else:
                      otu_filter_file.write('\n')
              
      
      
