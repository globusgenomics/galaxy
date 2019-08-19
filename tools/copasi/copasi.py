
from COPASI import *
import sys

def main(args):
  assert CCopasiRootContainer.getRoot() != None
  # create a datamodel
  dataModel = CCopasiRootContainer.addDatamodel()
  assert CCopasiRootContainer.getDatamodelList().size() == 1
  # the only argument to the main routine should be the name of an SBML file
  if len(sys.argv) == 11:

      filename = sys.argv[1]
      output_file = sys.argv[2]
      num_steps = int(sys.argv[3])
      final_time = float(sys.argv[4])
      algorithm = sys.argv[5]
      rel_tol = float(sys.argv[6])
      abs_tol = float(sys.argv[7])
      max_steps = int(sys.argv[8])
      use_rand_seed = bool(sys.argv[9]=="1")
      seed_value = int(sys.argv[10])

      try:
          # load the model 
          dataModel.importSBML(filename)
      except:
          print >> sys.stderr,  "Error while importing the model from file named \"" + filename + "\"." 
          return 1
      model = dataModel.getModel()
      assert model != None
      # create a report with the correct filename and all the species against
      # time.
      reports = dataModel.getReportDefinitionList()
      # create a report definition object
      report = reports.createReportDefinition("Report", "Output for timecourse")
      # set the task type for the report definition to timecourse
      report.setTaskType(CCopasiTask.timeCourse)
      # we don't want a table
      report.setIsTable(False)
      # the entries in the output should be seperated by a ", "
      report.setSeparator(CCopasiReportSeparator(", "))

      # we need a handle to the header and the body
      # the header will display the ids of the metabolites and "time" for
      # the first column
      # the body will contain the actual timecourse data
      header = report.getHeaderAddr()
      body = report.getBodyAddr()
      
      body.push_back(CRegisteredObjectName(CCopasiObjectName(dataModel.getModel().getCN().getString() + ",Reference=Time").getString()))
      body.push_back(CRegisteredObjectName(report.getSeparator().getCN().getString()))
      header.push_back(CRegisteredObjectName(CCopasiStaticString("time").getCN().getString()))
      header.push_back(CRegisteredObjectName(report.getSeparator().getCN().getString()))

      iMax = model.getMetabolites().size()
      for i in range(0,iMax):
          metab = model.getMetabolite(i)
          assert metab != None
          # we don't want output for FIXED metabolites right now
          if (metab.getStatus() != CModelEntity.FIXED):
              # we want the concentration oin the output
              # alternatively, we could use "Reference=Amount" to get the
              # particle number
              body.push_back(CRegisteredObjectName(metab.getObject(CCopasiObjectName("Reference=Concentration")).getCN().getString()))
              # add the corresponding id to the header
              header.push_back(CRegisteredObjectName(CCopasiStaticString(metab.getSBMLId()).getCN().getString()))
              # after each entry, we need a seperator
              if(i!=iMax-1):
                body.push_back(CRegisteredObjectName(report.getSeparator().getCN().getString()))
                header.push_back(CRegisteredObjectName(report.getSeparator().getCN().getString()))


      # get the trajectory task object
      trajectoryTask = dataModel.getTask("Time-Course")
      # if there isn't one
      if (trajectoryTask == None):
          # create a one
          trajectoryTask = CTrajectoryTask()
          # add the time course task to the task list
          # this method makes sure the object is now owned by the list
          # and that SWIG does not delete it
          dataModel.getTaskList().addAndOwn(trajectoryTask)


      # run a deterministic time course
      if algorithm=="deterministic":
        trajectoryTask.setMethodType(CCopasiMethod.deterministic)
      elif algorithm=="stochastic":
        trajectoryTask.setMethodType(CCopasiMethod.stochastic)
      else:
        print >> sys.stderr, "Internal error: invalid algorithm." 
        return 1

      # pass a pointer of the model to the problem
      trajectoryTask.getProblem().setModel(dataModel.getModel())

      # actiavate the task so that it will be run when the model is saved
      # and passed to CopasiSE
      trajectoryTask.setScheduled(True)

      # set the report for the task
      trajectoryTask.getReport().setReportDefinition(report)
      # set the output filename
      trajectoryTask.getReport().setTarget(output_file)
      # don't append output if the file exists, but overwrite the file
      trajectoryTask.getReport().setAppend(False)

      # get the problem for the task to set some parameters
      problem = trajectoryTask.getProblem()

      # simulate num_steps steps
      problem.setStepNumber(num_steps)
      # start at time 0
      dataModel.getModel().setInitialTime(0.0)
      # simulate a duration of final_time time units
      problem.setDuration(final_time)
      # tell the problem to actually generate time series data
      problem.setTimeSeriesRequested(True)

      # set some parameters for the LSODA method through the method
      method = trajectoryTask.getMethod()


      if algorithm=="deterministic":
        parameter = method.getParameter("Relative Tolerance")
        assert parameter != None
        assert parameter.getType() == CCopasiParameter.UDOUBLE
        parameter.setValue(rel_tol)
        parameter = method.getParameter("Absolute Tolerance")
        assert parameter != None
        assert parameter.getType() == CCopasiParameter.UDOUBLE
        parameter.setValue(abs_tol)
        parameter = method.getParameter("Max Internal Steps")
        assert parameter != None
        assert parameter.getType() == CCopasiParameter.UINT
        parameter.setValue(max_steps)
      elif algorithm=="stochastic":
        parameter = method.getParameter("Max Internal Steps")
        assert parameter != None
        assert parameter.getType() == CCopasiParameter.INT
        parameter.setValue(max_steps)
        parameter = method.getParameter("Use Random Seed")
        assert parameter != None
        assert parameter.getType() == CCopasiParameter.BOOL
        parameter.setValue(use_rand_seed)
        parameter = method.getParameter("Random Seed")
        assert parameter != None
        assert parameter.getType() == CCopasiParameter.UINT
        parameter.setValue(seed_value)
        

      result=True
      try:
          # now we run the actual trajectory
          result=trajectoryTask.process(True)
      except:
          # check if there are additional error messages
          if CCopasiMessage.size() > 0:
              # print the messages in chronological order
              print >> sys.stderr, CCopasiMessage.getAllMessageText(True)
          print >> sys.stderr,  "Error. Running the time course simulation failed." 
          return 1
      if result == False:
          # check if there are additional error messages
          if CCopasiMessage.size() > 0:
              # print the messages in chronological order
              print >> sys.stderr, CCopasiMessage.getAllMessageText(True)
          print >> sys.stderr,  "Error. Running the time course simulation failed." 
          return 1

      # look at the timeseries
      timeSeries = trajectoryTask.getTimeSeries()
      # we simulated num_steps steps, including the initial state, this should be
      # num_steps+1 steps in the timeseries
      assert timeSeries.getRecordedSteps() == num_steps+1
      print  "The time series consists of " , timeSeries.getRecordedSteps() , "." 
      print  "Each step contains " , timeSeries.getNumVariables() , " variables." 
      print  "The final state is: " 
      iMax = timeSeries.getNumVariables()
      lastIndex = timeSeries.getRecordedSteps() - 1
      for i in range(0,iMax):
          # here we get the particle number (at least for the species)
          # the unit of the other variables may not be particle numbers
          # the concentration data can be acquired with getConcentrationData
          print timeSeries.getTitle(i) + ": " , timeSeries.getData(lastIndex, i) 
      # the CTimeSeries class now has some new methods to get all variable titles
      # as a python list (getTitles())
      # and methods to get the complete time course data for a certain variable based on
      # the variables index or the corresponding model object.
      # E.g. to get the particle numbers of the second variable as a python list
      # you can use getDataForIndex(1) and to get the concentration data you use
      # getConcentrationDataForIndex(1)
      # To get the complete particle number data for the second metabolite of the model
      # you can use getDataForObject(model.getMetabolite(1)) and to get the concentration
      # data you use getConcentrationDataForObject.
      #print timeSeries.getTitles()
      #print timeSeries.getDataForIndex(1)
      #print timeSeries.getDataForObject(model)

      #replace default report delimiter (", ") with tab
      s=open(output_file).read()
      s=s.replace(", ", "\t")
      f=open(output_file, 'w')
      f.write(s)
      f.flush()
      f.close()

  else:
      print >> sys.stderr, "Internal error: Usage: copasi.py $sbml_input $out_file $num_steps $duration $algorithm $options.rel_tol $options.abs_tol $options.max_steps $options.stoch_options.use_rand_seed $options.stoch_options.seed_value" 
      return 1;

if(__name__ == '__main__'):
   main(sys.argv[1:]) 
