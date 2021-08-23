#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2021.1.2),
    on Mon Mar 22 11:15:44 2021
    and then heavily modified by Alex, Liam and Josh
"""
from __future__ import absolute_import, division
from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors
from psychopy import event
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np  # whole numpy lib is available, prepend 'np.'
import pylab #for some frametimes plotting
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os, shutil  # handy system and path functions
import sys  # to get file system encoding
import pandas as pd

from psychopy.hardware import keyboard

from datetime import datetime

debug=False #Print more information to console
autopilot=False
demo=False

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

#Set up screen stuff
frameTolerance = 0.001  # how close to onset before 'same' frame
refreshRatePlanningOn = 60

monitorwidth = 28.5 #monitor width in cm
viewdist = 57 #cm
scrn=0 #0 to use main screen (or second screen with mirroring), 1 to use external screen connected to computer
widthPix = 2880
heightPix = 1800
bgColor = [0.004,0.004,0.004]
fullscr = 1 #Full screen necessary for good timing

#monitor characteristcs
import psychopy
from psychopy import monitors
monitorname = 'testmonitor'
waitBlank = False
mon = monitors.Monitor(monitorname,width=monitorwidth, distance=viewdist)#relying on  monitorwidth cm (39 for Mitsubishi to do deg calculations) and gamma info in calibratn
mon.setSizePix( (widthPix,heightPix) )

def openMyStimWindow(): #make it a function because have to do it several times, want to be sure is identical each time
    myWin = visual.Window(monitor=mon,size=(widthPix,heightPix),allowGUI=False,units='pix',color=bgColor,colorSpace='rgb',fullscr=fullscr,
              screen=scrn,waitBlanking=waitBlank, useFBO=True, winType='pyglet' ) #pygame doesn't work, don't know why. Works in textLocationTest.py
    return myWin
myWin = openMyStimWindow()

from psychopy import info #20 June 2021: psychopy.info doesn't work even though psychopy.visual does https://discourse.psychopy.org/t/psychopy-info-disappeared/23619/2
checkRefreshEtc = True
refreshMsg2 = ''
if not checkRefreshEtc:
    refreshMsg1 = 'REFRESH RATE WAS NOT CHECKED'
    refreshRateWrong = False
else: #checkRefreshEtc
    runInfo = info.RunTimeInfo(
            win=myWin,    ## a psychopy.visual.Window() instance; None = default temp window used; False = no win, no win.flips()
            refreshTest='grating', ## None, True, or 'grating' (eye-candy to avoid a blank screen)
            verbose=True, ## True means report on everything 
            userProcsDetailed=True  ## if verbose and userProcsDetailed, return (command, process-ID) of the user's processes
            )
    #print(runInfo)
    #logging.info(runInfo)
    #check screen refresh is what assuming it is ##############################################
    Hzs=list()
    myWin.flip(); myWin.flip();myWin.flip();myWin.flip();
    myWin.setRecordFrameIntervals(True) #otherwise myWin.fps won't work
    for i in range(50):
        myWin.flip()
        Hzs.append( myWin.fps() )  #varies wildly on successive runs!
    myWin.setRecordFrameIntervals(False)
    # end testing of screen refresh########################################################
    Hzs = np.array( Hzs );     Hz= np.median(Hzs)
    msPerFrame= 1000./Hz
    refreshMsg1= 'Frames per second ~='+ str( np.round(Hz,1) )
    refreshRateTolerancePct = 3
    refreshRateObserved = np.median(Hzs)
    pctOff = abs( (refreshRateObserved-refreshRatePlanningOn) / refreshRatePlanningOn)
    refreshRateWrong =  pctOff > (refreshRateTolerancePct/100.)
    if refreshRateWrong:
        refreshMsg1 += ' BUT'
        refreshMsg1 += ' program assumes ' + str(refreshRatePlanningOn)
        refreshMsg2 =  'which is off by more than' + str(round(refreshRateTolerancePct,0)) + '%!!'
    else:
        refreshMsg1 += ', which is close enough to desired val of ' + str( round(refreshRatePlanningOn,1) )
    myWinRes = myWin.size
    myWin.allowGUI =True
myWin.close() #have to close window to show dialog box

# Store info about the experiment session
from psychopy import __version__ # Get the PsychoPy version currently in use
psychopyVersion = __version__
expName = 'noiseMot_exp2_noise'  # from the Builder filename that created this script
expInfo = {'motbox_path': '../', 'participant': '999'}

myDlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName, show=False)

myDlg.addText(refreshMsg1, color='Black')
if refreshRateWrong:
    myDlg.addText(refreshMsg2, color='Red')
dimGreyForDlgBox = 'DimGrey'

if checkRefreshEtc and (myWinRes != [widthPix,heightPix]).any():
    msgWrongResolution = 'Screen apparently NOT the desired resolution of '+ str(widthPix)+'x'+str(heightPix)+ ' pixels!'+ \
                            '\ninstead it is ' + str(myWinRes)
    myDlg.addText(msgWrongResolution, color='Red')
    logging.error(msgWrongResolution)
myDlg.addText('Note: to abort press ESC at a trials response screen', color=dimGreyForDlgBox) 

flaggedProcessesMsg = False
if runInfo['systemUserProcFlagged'] is not None:
    if len(runInfo['systemUserProcFlagged']):
        flaggedProcessesMsg = 'Other programs running: (command, process-ID)' + str(runInfo['systemUserProcFlagged'])
if flaggedProcessesMsg:
    myDlg.addText(flaggedProcessesMsg, color='Red')
myDlg.show()    

if myDlg.OK == False:
    core.quit()  # user pressed cancel
if refreshRateWrong:
    logging.error(refreshMsg1+refreshMsg2)
else: logging.info(refreshMsg1+refreshMsg2)

expInfo['expName'] = expName
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
if not demo and not debug: #Save a copy of the code so know exactly what the code was for each participant
    shutil.copy2(sys.argv[0], filename + '.py') #https://stackoverflow.com/questions/123198/how-can-a-file-be-copied

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=runInfo,
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.DATA)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
#End screen and window testing
import matplotlib
def plotFrameIntervals(intervals_msec):
        matplotlib.use('Qt5Agg')  # change this to control the plotting 'back end'
        m = pylab.mean(intervals_msec)
        sd = pylab.std(intervals_msec)        
        msg = "Mean=%.1fms, s.d.=%.2f, 99%%CI(frame)=%.2f-%.2f"
        distString = msg % (m, sd, m - 2.58 * sd, m + 2.58 * sd)
        nTotal = len(intervals_msec)
        nDropped = sum(intervals_msec > longFrameLimit)
        msg = "Across all trials, Dropped / Frames = %i/%i = %.3f%%"
        droppedString = msg % (nDropped, nTotal, 100 * nDropped / float(nTotal))
        # plot the frameintervals
        pylab.figure(figsize=[12, 8])
        pylab.subplot(1, 2, 1)
        pylab.plot(intervals_msec, '-')
        pylab.ylabel('t (ms)')
        pylab.xlabel('frame N')
        pylab.title(droppedString)
        
        pylab.subplot(1, 2, 2)
        pylab.hist(intervals_msec, 50, histtype='stepfilled')
        pylab.xlabel('t (ms)')
        pylab.ylabel('n frames')
        pylab.title(distString)
        pylab.show()
# setup the window for the actual practice trials
win = openMyStimWindow()
win.setRecordFrameIntervals(False) 
frameTimeTolerance=.2 #proportion longer than refreshRate that will not count as a miss
# Ay refresh that takes longer than refreshThreshold will be considered a "dropped"
# frame and increase the count of win.nDroppedFrames, during periods when win.recordFrameIntervals = True
win.refreshThreshold = 1/refreshRateObserved + 0.20*(1/refreshRateObserved)
# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "welcome"
welcomeClock = core.Clock()
text_2 = visual.TextStim(win=win, name='text_2',
    text='In each trial, there will be 2 objects, which will move towards the center of the screen, then disappear. You will then need to press DOWN if you thought they disappeared before passing the red cross, or UP if you thought the objects disappeared beyond the red cross.\n\nKeep your eyes on the red cross whenever it is visible.\n\nPress space to continue',
    font='Arial',
    pos=(0, 0), height=20, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
key_resp_3 = keyboard.Keyboard()

# Initialize components for Routine "intro_block"
intro_blockClock = core.Clock()
text_3 = visual.TextStim(win=win, name='text_3',
    text='First, there will be 8 practice trials.\n\nPress space to continue',
    font='Arial',
    pos=(0, 0), height=20, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

key_resp_4 = keyboard.Keyboard()
# Initialize components for Routine "trial"
trialClock = core.Clock()
sys.path.insert(0, expInfo['motbox_path']) # to load module from neighbouring folder
import motbox
import time
import copy

nTargets = 2 # number of targets (that will be highlighted in cue phase)
n_objects = 2 # number of objects (could be loaded automatically from file

# configure experiment
cue_time = 2     # duration of cueing phase in seconds

# hide mouse cursor
win.mouseVisible = False

# Initialize components for Routine "before_trials"
before_trialsClock = core.Clock()
text_5 = visual.TextStim(win=win, name='text_5',
    text='Signal to the experimenter when you are ready to begin the real trials.',
    font='Arial',
    pos=(0, 0), height=20, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
key_resp_2 = keyboard.Keyboard()
# Initialize components for Routine "fix"
fixClock = core.Clock()
text = visual.TextStim(win=win, name='text',
    text='+',
    font='Arial',
    pos=(0, 5), height=100, wrapWidth=None, ori=0, 
    color='red', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-10);

# Initialize components for Routine "trial"
trialClock = core.Clock()


sys.path.insert(0, expInfo['motbox_path']) # to load module from neighbouring folder
import motbox
import time
import copy

# configure experiment
#trial_length = 6 # duration of trial in seconds
#cue_time = 2     # duration of cueing phase in seconds LJM trial_length and cue_time are defined later

# hide mouse cursor
win.mouseVisible = False

# Initialize components for Routine "bye"
byeClock = core.Clock()
text_4 = visual.TextStim(win=win, name='text_4',
    text='Thanks!',
    font='Arial',
    pos=(0, 0), height=20, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

noise = visual.NoiseStim(
    win=win, name='noise',
    mask=None,
    ori=0.0, pos=(0, 0), size=(2048,2048),
    sf=None,
    color=[1,1,1], colorSpace='rgb', opacity=None, blendmode='avg', contrast=1.0,
    texRes=128, filter=None,
    noiseType='Uniform', noiseElementSize= (8, 8),
    interpolate=False, depth=0.0, units = 'pix')
noise.buildNoise()

# ------Prepare to start Routine "welcome"-------
continueRoutine = True
# update component parameters for each repeat
key_resp_3.keys = []
key_resp_3.rt = []
_key_resp_3_allKeys = []
# keep track of which components have finished
welcomeComponents = [text_2, key_resp_3]
for thisComponent in welcomeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
welcomeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "welcome"-------
while continueRoutine:
    # get current time
    t = welcomeClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=welcomeClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_2* updates
    if text_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_2.frameNStart = frameN  # exact frame index
        text_2.tStart = t  # local t and not account for scr refresh
        text_2.tStartRefresh = tThisFlipGlobal  # on global time
        text_2.setAutoDraw(True)
    
    # *key_resp_3* updates
    waitOnFlip = False
    if key_resp_3.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp_3.frameNStart = frameN  # exact frame index
        key_resp_3.tStart = t  # local t and not account for scr refresh
        key_resp_3.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp_3, 'tStartRefresh')  # time at next scr refresh
        key_resp_3.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp_3.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp_3.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp_3.status == STARTED and not waitOnFlip:
        theseKeys = key_resp_3.getKeys(keyList=['space'], waitRelease=False)
        _key_resp_3_allKeys.extend(theseKeys)
        if len(_key_resp_3_allKeys):
            key_resp_3.keys = _key_resp_3_allKeys[-1].name  # just the last key pressed
            key_resp_3.rt = _key_resp_3_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in welcomeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "welcome"-------
for thisComponent in welcomeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# the Routine "welcome" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

practiceConditionInfo={}
practiceConditionInfo['conditions'] = ['500x500xdynamic','500x1000xdynamic','1000x500xdynamic',
    '1000x1000xdynamic','500x500xstatic','500x1000xstatic','1000x500xstatic','1000x1000xstatic']
practiceConditionInfo['nTrials'] = 1

pracStairs=[]
for thisCondition in practiceConditionInfo['conditions']:
    thisStair = data.StairHandler(startVal=0, 
        stepSizes=[4,4,2,2,1,1], 
        nTrials=5, nUp=1, nDown=1, 
        extraInfo=thisCondition, 
        stepType='lin', minVal=-4, 
        maxVal=4)
    pracStairs.append(thisStair)

# ------Prepare to start Routine "intro_block"-------
continueRoutine = True
# update component parameters for each repeat
win.blendMode = 'add'
key_resp_4.keys = []
key_resp_4.rt = []
_key_resp_4_allKeys = []
# keep track of which components have finished
intro_blockComponents = [text_3, key_resp_4]
for thisComponent in intro_blockComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
intro_blockClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "intro_block"-------
while continueRoutine:
    # get current time
    t = intro_blockClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=intro_blockClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_3* updates
    if text_3.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_3.frameNStart = frameN  # exact frame index
        text_3.tStart = t  # local t and not account for scr refresh
        text_3.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_3, 'tStartRefresh')  # time at next scr refresh
        text_3.setAutoDraw(True)
    
    # *key_resp_4* updates
    waitOnFlip = False
    if key_resp_4.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp_4.frameNStart = frameN  # exact frame index
        key_resp_4.tStart = t  # local t and not account for scr refresh
        key_resp_4.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp_4, 'tStartRefresh')  # time at next scr refresh
        key_resp_4.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp_4.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp_4.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp_4.status == STARTED and not waitOnFlip:
        theseKeys = key_resp_4.getKeys(keyList=['space'], waitRelease=False)
        _key_resp_4_allKeys.extend(theseKeys)
        if len(_key_resp_4_allKeys):
            key_resp_4.keys = _key_resp_4_allKeys[-1].name  # just the last key pressed
            key_resp_4.rt = _key_resp_4_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in intro_blockComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "intro_block"-------
for thisComponent in intro_blockComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# the Routine "intro_block" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
practice_trials = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=[],
    seed=None, name='practice_trials')
thisExp.addLoop(practice_trials)  # add the loop to the experiment
thisPractice_trial = practice_trials.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisPractice_trial.rgb)
if thisPractice_trial != None:
    for paramName in thisPractice_trial:
        exec('{} = thisPractice_trial[paramName]'.format(paramName))

# This is a template for one object
o1 = visual.Rect(
    win=win, name ='o1',
    units="pix",
    size=(50,100),
    fillColor=['red'],
    lineColor=['red'], depth =-1.0)
o1.setAutoDraw(False) #Because copies will be made of it by Puppetteer (P), and it will draw them, so o1 is never drawn

for trialN in range(practiceConditionInfo['nTrials']):
    shuffle(pracStairs)
    for thisStair in pracStairs:
        thisIntensity = next(thisStair)
        print(thisIntensity)
        track_filename = 'trajectories/T%00d%s.csv' % (thisIntensity,thisStair.extraInfo)
        if ((trialN + 1) % 2) > 0:
            track_filename = 'trajectories/T%00d%sf.csv' % (thisIntensity,thisStair.extraInfo)
        currentLoop = practice_trials
        # abbreviate parameter names if possible (e.g. rgb = thisPractice_trial.rgb)
        if thisPractice_trial != None:
            for paramName in thisPractice_trial:
                exec('{} = thisPractice_trial[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial"-------
        continueRoutine = True
        
        # initialize motbox trajectory
        T = motbox.Track()
        T.load_from_csv(track_filename, delim = ',') # load it from file
        
        # scale it with given pixels per inch. This needs to be set for proper values
        T.scale(1) # pixels per unit in the trajectory file. Gets multiplied by the trajectory files 
        
        # Puppeteer handles the actual moving of objects
        P = motbox.Puppeteer()
        P.track = T
        # this creates multiple copies from the template
        P.clone_template_psychopy(o1, n_objects)
        
        # Initialize object state as participant not having clicked on it
        for obj in P.objects:
            obj.pressed   = False
        
        win.blendMode = 'avg' # 'add' #set proper blending mode JJC Average (avg) places the new stimulus over the old one with a transparency given by its opacity.
        iframe = 1 # this is iterator through frames
        
        size_blink  = 8       # how many frames are there between size change 
        size_big_coef    = 2  # how large are objects when their size is increased
        size_normal_coef = 1  # normal size coeficient, should be 1
        size_normal = P.objects[0].size
        size_coef = size_normal_coef 
        
        # set start time
        # in this particular study, we are not randomizing starts of the trials
        # to randomize starts, set variable start_time2 to value given value
        nRowsForTrajectoryFile = T.x.shape[0] #number of rows in trajectory file gives the length of this trial
        timeGrain = T.time[1] #second cell of time column gives the time grain for this trial
        trial_length = (nRowsForTrajectoryFile - 1) * timeGrain 
        
        start_time2 = 0
        stop_time = cue_time + trial_length 
        
        P.update_positions_psychopy(start_time2)
        pos_obj1 = P.objects[0].pos
        
        # keep track of which components have finished
        trialComponents = [P.objects[0], noise]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        P.objects[0].finalx0 = -999 #this is just for recording the final position to the data file
        P.objects[0].finaly0 = -999
        P.objects[0].penultimatex0 = -999 #this is just for recording the penultimate position to the data file
        P.objects[0].penultimatey0 = -999 #this is just for recording the penultimate position to the data file
        P.objects[0].antepenultimatex0 = -999 #this is just for recording the antepenultimate position to the data file
        P.objects[0].antepenultimatey0 = -999 #this is just for recording the antepenultimate position to the data file
        P.objects[0].preantepenultimatex0 = -999 #this is just for recording the preantepenultimate position to the data file
        P.objects[0].preantepenultimatey0 = -999 #this is just for recording the preantepenultimate position to the data file
        time.sleep(0.2)
        # reset timers
        t = 0
        nKeys = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        finishedCriticalStimuli = False
        
        # -------Run Routine "trial"-------
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=trialClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            #Calculate offset between global clock and time since stimuli started,
            #because what gets recorded by win.timeOnFlip is from global clock, not trial clock
            #print('time relative to stimulus start =',t, 'time according to tThisFlipGlobal=',tThisFlipGlobal)
            #Problem is that tThisFlipGlobal is expected time of next flip,but I need current time
            #print('trial starting at approximately tThisFlipGlobal:', tThisFlipGlobal, 'tThisFlip (from trialClock):',tThisFlip)
            #From the psychopy code in getFutureFlipTime, it seems like logging.defaultClock is the global clock and it
            #calculates the difference using logging.defaultClock.getLastResetTime() - clock.getLastResetTime()
            timeSuspectGlobalClockAheadBy = tThisFlipGlobal - t
            globalClockAheadBy = trialClock.getLastResetTime() - logging.defaultClock.getLastResetTime()
            #print('Based on win.getFutureFlipTime, Expecting global clock to be ahead by a bit less than',timeSuspectGlobalClockAheadBy,
            #        'and actually seems globalClockAheadBy=',globalClockAheadBy, 'so will subtract that from the time recorded to file')
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            
            noise.draw( )
            o1.setAutoDraw(False)
            #cue.setAutoDraw(False)
            # cueing of the targets phase at beginning of trial
            if t >= 0 and t <= cue_time:
                text.draw( )
                # simple iterator through frames
                # we are not using t, as it is not an integer (we could round though)
                iframe = iframe + 1     
                
                # update and draw objects
                P.update_positions_psychopy(start_time2)    
                for i2 in range(n_objects):
                    P.objects[i2].draw() 
                
            # move phase
            if t >= cue_time and t <= stop_time:
                event.clearEvents()
                text.draw( )
                # normalize sizes 
                for i2 in range(nTargets):
                    P.objects[i2].size = size_normal
                # update and draw objects
                P.update_positions_psychopy(start_time2 + t-cue_time)
                for i2 in range(n_objects):
                    P.objects[i2].draw() 
                numOfFinalFramesToRecord = 40
                tLastFramesStart = stop_time - numOfFinalFramesToRecord*(1.0/refreshRateObserved) #time after which should be the last few frames
                if t > tLastFramesStart: #make sure pick up the last frames
                    #Record the time of the very last frame that the stimuli are still on, by doing it over and over, so when the code ceases it will be correct.
                    P.objects[0].tLastFrame = -999 #dummy value which will hopefully get overwritten on win flip by the next line, over and over until the very last frame
                    win.timeOnFlip(P.objects[0], 'tLastFrame')  # set P.objects[0].tLastFrame to time at next scr refresh
                    
                     #To record the penultimate position to the data file, on each pass-through, set it to final before setting final to current, so that
                    #at the last frame, final will be final
                    P.objects[0].preantepenultimatex0 = P.objects[0].antepenultimatex0
                    P.objects[0].preantepenultimatey0 = P.objects[0].antepenultimatey0
                    P.objects[0].antepenultimatex0 = P.objects[0].penultimatex0
                    P.objects[0].antepenultimatey0 = P.objects[0].penultimatey0
                    P.objects[0].penultimatex0 = P.objects[0].finalx0 #You might think could rely on whatever the final x is, but this making sure recording the final x *drawn*
                    P.objects[0].penultimatey0 = P.objects[0].finaly0 #You might think could rely on whatever the final x is, but this making sure recording the final x *drawn*
                    
                    #Record last x,y of object 0
                    P.objects[0].finalx0 = P.objects[0].pos[0] #You might think could rely on whatever the final x is, but this making sure recording the final x *drawn*
                    P.objects[0].finaly0 = P.objects[0].pos[1]
    
                    #But win.timeOnFlip won't work for recording a whole list or array of times because timeOnFlip can't append to an array.
                    #So for recording the interframe intervals of these times the conventional way is https://www.psychopy.org/general/timing/detectingFrameDrops.html
                    win.recordFrameIntervals = True
                    
                    #Another way is to record it directly myself
                    #https://github.com/psychopy/psychopy/blob/e45520f446697d5b5fad035fa0e4e50088ffa535/psychopy/visual/window.py
                    #That means calling my own function that gets the time of the flip.
                    #To get the time of the flip, it seems that the Psychopy functions use win._frameTime, which is set by logging.defaultClock.getTime() automatically
                    #Alternatively, I could probably use the automatic log of win._frameTimes, but then at the time of calling addData, to know 
                    #which are the critical frames, I'll need to record how long the array of win._frameTimes is at the last frame.
                    #win.callOnFlip(appendToLastFrameTimes, frameTime=win._frameTime
                    #win.timeOnFlip( lastFrameTimes[0], 
    
            # query phase
            if t > stop_time:
                text.draw( )
                keys = event.getKeys()
                for thisKey in keys:
                    if (thisKey == 'up'):
                        wasUp = 1
                        nKeys = nKeys + 1
                    elif (thisKey == 'down'): 
                        wasUp = 0
                        nKeys = nKeys + 1
                if not finishedCriticalStimuli:
                    lenFrameTimesUntilStimuliFinish = len(win._frameTimes)
                finishedCriticalStimuli = True
                win.recordFrameIntervals = False
                if 'dynamic' in thisStair.extraInfo:
                    if noise.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                        noise.frameNStart = frameN  # exact frame index
                        noise.tStart = t  # local t and not account for scr refresh
                        noise.tStartRefresh = tThisFlipGlobal  # on global time, but not as accurate as win.time so overrule it in next line with that
                        win.timeOnFlip(noise, 'tStartRefresh')  # time at next scr refresh
                        noise.status = STARTED #AOH
    
                    if noise.status == STARTED:
                        if noise._needBuild:
                            noise.buildNoise()
                        else:
                            noise.draw() #AOH
                            if (frameN-noise.frameNStart) % 2 == 0:
                                noise.updateNoise()
                else:
                    noise.draw( )
                if nKeys > 0:
                    thisStair.addResponse(wasUp)
                    expInfo['response'] = keys
                    expInfo['wasUp'] = wasUp
                    expInfo['condition'] = thisStair.extraInfo
                    if 'f' in track_filename:
                        expInfo['flipStatus'] = 'flipped'
                    else: expInfo['flipStatus'] = 'not_flipped'
                    expInfo['offset'] = -thisIntensity
                    expInfo['trialNumber'] = 'practice'
                    continueRoutine = False
                    for i2 in range(n_objects):
                        P.objects[i2].setAutoDraw(False)
                    thisExp.saveAsWideText(filename+'.csv', fileCollisionMethod='overwrite',delim='auto')
                    time.sleep(0.05) 
            
            # *o1* updates
            if P.objects[0].status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                P.objects[0].frameNStart = frameN  # exact frame index
                P.objects[0].tStart = t  # local t and not account for scr refresh
                win.timeOnFlip(P.objects[0], 'tStartRefresh')  # time at next scr refresh
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # save mouse data
        win.blendmode = 'avg'
        if debug:
            print('globalClockAheadBy=',globalClockAheadBy)
            print('P.objects[0].tStartRefresh - globalClockAheadBy = ',P.objects[0].tStartRefresh - globalClockAheadBy)
        practice_trials.addData('P.objects[0].started', P.objects[0].tStartRefresh - globalClockAheadBy)
        practice_trials.addData('P.objects[0].tLastFrame', P.objects[0].tLastFrame - globalClockAheadBy)
        practice_trials.addData('P.objects[0].penultimatex0',P.objects[0].penultimatex0)
        practice_trials.addData('P.objects[0].penultimatey0',P.objects[0].penultimatey0)
        practice_trials.addData('P.objects[0].antepenultimatex0',P.objects[0].antepenultimatex0)
        practice_trials.addData('P.objects[0].antepenultimatey0',P.objects[0].antepenultimatey0)
        practice_trials.addData('P.objects[0].preantepenultimatex0',P.objects[0].preantepenultimatex0)
        practice_trials.addData('P.objects[0].preantepenultimatey0',P.objects[0].preantepenultimatey0)
        practice_trials.addData('P.objects[0].finalx0',P.objects[0].finalx0)
        practice_trials.addData('P.objects[0].finaly0',P.objects[0].finaly0)
        
    #Record final frametimes when critical stimuli on
        practice_trials.addData('win.nDroppedFrames',win.nDroppedFrames)
        intervals_msec = pylab.array(win.frameIntervals) * 1000
        intervals_msec_last_thistrial = intervals_msec[-numOfFinalFramesToRecord:] #last numOfFinalFramesToRecord
        #intervals_msec = intervals_msec[numOfFinalFramesToRecord
        if debug:
            print('Frame times up to stimuli cessation according to win.frameIntervals=', intervals_msec_last_thistrial)
            print('In last',numOfFinalFramesToRecord,' frames across all trials ', win.nDroppedFrames, ' frames were dropped according to nDroppedFrames.')
        practice_trials.addData('finalStimFrameTimes',intervals_msec_last_thistrial)
        worstLastFramesDur = max(intervals_msec_last_thistrial)
        practice_trials.addData('worstLastFramesDur',worstLastFramesDur)
        #Count the number of timing hiccups and print out and save to log some information about them
        longFrameLimit = np.round(1000/refreshRateObserved*(1.0+frameTimeTolerance),2)
        idxsInterframeLong = np.where( np.array(intervals_msec_last_thistrial) > longFrameLimit ) [0] #frames that exceeded frameTimeTolerance of expected duration
        numCasesInterframeLong = len( idxsInterframeLong )
        practice_trials.addData('timingHiccupsInLastFramesOfStimuli',numCasesInterframeLong)
        if numCasesInterframeLong >0 and (not demo):
           longFramesStr =  'ERROR,'+str(numCasesInterframeLong)+' frames were longer than '+str(longFrameLimit)+' ms'
           if demo: 
             longFramesStr += 'not printing them all because in demo mode'
           else:
               longFramesStr += ' apparently screen refreshes skipped, interframe durs were:'+\
                        str( np.around(  intervals_msec_last_thistrial[idxsInterframeLong] ,1  ) )+ ' and was these frames: '+ str(idxsInterframeLong)
           if longFramesStr != None:
                    logging.error( 'trialnum='+str(practice_trials.thisN)+' '+longFramesStr )
                    if not demo:
                        flankingAlso=list()
                        for idx in idxsInterframeLong: #also print timing of one before and one after long frame
                            if idx-1>=0:
                                flankingAlso.append(idx-1)
                            else: flankingAlso.append(np.NaN)
                            flankingAlso.append(idx)
                            if idx+1<len(intervals_msec_last_thistrial):  flankingAlso.append(idx+1)
                            else: flankingAlso.append(np.NaN)
                        flankingAlso = np.array(flankingAlso)
                        flankingAlso = flankingAlso[np.logical_not(np.isnan(flankingAlso))]  #remove nan values
                        flankingAlso = flankingAlso.astype(np.integer) #cast as integers, so can use as subscripts
                        logging.info( 'flankers also='+str( np.around( intervals_msec_last_thistrial[flankingAlso], 1) )  ) #because this is not an essential error message, as previous one already indicates error
                        #As INFO,it won't fill up the console when console set to WARNING or higher
    
        if debug and (practice_trials.thisN==2): #Quit and Plot frame intervals
            win.close()
            plotFrameIntervals(intervals_msec)
            core.quit()
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
    
# completed 1 repeats of 'practice_trials'

#REAL TRIALS

conditionInfo={}
conditionInfo['conditions'] = ['500x500xdynamic','500x1000xdynamic','1000x500xdynamic',
    '1000x1000xdynamic','500x500xstatic','500x1000xstatic','1000x500xstatic','1000x1000xstatic']
conditionInfo['nTrials'] = 125

stairs=[]
for thisCondition in conditionInfo['conditions']:
    thisStartVal = randint(-4,4)
    thisStair = data.StairHandler(startVal=thisStartVal, 
        stepSizes=[4,4,2,2,1,1], 
        nTrials=125, nUp=1, nDown=1, 
        extraInfo=thisCondition, 
        stepType='lin', minVal=-4, 
        maxVal=4)
    stairs.append(thisStair)
    

# ------Prepare to start Routine "before_trials"-------
win.blendMode = 'add'
cue_time = 0.5
continueRoutine = True
# update component parameters for each repeat
key_resp_2.keys = []
key_resp_2.rt = []
_key_resp_2_allKeys = []
# keep track of which components have finished
before_trialsComponents = [text_5, key_resp_2]
for thisComponent in before_trialsComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
before_trialsClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "before_trials"-------
while continueRoutine:
    # get current time
    t = before_trialsClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=before_trialsClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_5* updates
    if text_5.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_5.frameNStart = frameN  # exact frame index
        text_5.tStart = t  # local t and not account for scr refresh
        text_5.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_5, 'tStartRefresh')  # time at next scr refresh
        text_5.setAutoDraw(True)
    
    # *key_resp_2* updates
    waitOnFlip = False
    if key_resp_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp_2.frameNStart = frameN  # exact frame index
        key_resp_2.tStart = t  # local t and not account for scr refresh
        key_resp_2.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp_2, 'tStartRefresh')  # time at next scr refresh
        key_resp_2.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp_2.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp_2.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp_2.status == STARTED and not waitOnFlip:
        theseKeys = key_resp_2.getKeys(keyList=['y'], waitRelease=False)
        _key_resp_2_allKeys.extend(theseKeys)
        if len(_key_resp_2_allKeys):
            key_resp_2.keys = _key_resp_2_allKeys[-1].name  # just the last key pressed
            key_resp_2.rt = _key_resp_2_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in before_trialsComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "before_trials"-------
for thisComponent in before_trialsComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# the Routine "before_trials" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
trials = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=[],
    seed=None, name='trials')
thisExp.addLoop(trials)  # add the loop to the experiment
thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
if thisTrial != None:
    for paramName in thisTrial:
        exec('{} = thisTrial[paramName]'.format(paramName))

o1.setAutoDraw(False) #Because copies will be made of it by Puppetteer (P), and it will draw them, so o1 is never drawn
# and circular highlighting
# this is related to current experiment (see poster)
cue = visual.Rect(
    win=win, name='cue',units='norm', 
    size=(1, 1),
    ori=0, pos=(0.5, 0.5),
    lineWidth=0, lineColor=[-1,1,-1], lineColorSpace='rgb',
    fillColor=None, fillColorSpace='rgb', autoDraw=False,
    opacity=1, depth=-1.0, interpolate=True)

trialNumber = 0

for trialN in range(conditionInfo['nTrials']):
    shuffle(stairs)
    for thisStair in stairs:
        thisIntensity = next(thisStair)
        track_filename = 'trajectories/T%00d%s.csv' % (thisIntensity,thisStair.extraInfo)
        if ((trialN + 1) % 2) > 0:
            track_filename = 'trajectories/T%00d%sf.csv' % (thisIntensity,thisStair.extraInfo)
        currentLoop = trials
        frameN = -1
        win.blendMode = 'avg'
     
        # ------Prepare to start Routine "trial"-------
        continueRoutine = True
        win.mouseVisible = False
        
        # initialize motbox trajectory
        T = motbox.Track()
        T.load_from_csv(track_filename, delim = ',') # load it from file
        
        # scale it with given pixels per inch. This needs to be set for proper values
        T.scale(1) # pixels per trajectory file coordinate unit 
        
        # Puppeteer handles the actual moving of objects
        P = motbox.Puppeteer()
        P.track = T
        # this creates multiple copies from the template
        P.clone_template_psychopy(o1, n_objects)
        
        # Initialize object state as participant not having clicked on it
        for obj in P.objects:
            obj.pressed   = False
        
        win.blendMode = 'avg' # set proper blending mode
        iframe = 1 # this is iterator through frames
        
        size_blink  = 8       # how many frames are there between size change 
        size_big_coef    = 2  # how large are objects when their size is increased
        size_normal_coef = 1  # normal size coeficient, should be 1
        size_normal = P.objects[0].size
        size_coef = size_normal_coef 
        
        # set start time
        # in this particular study, we are not randomizing starts of the trials
        # to randomize starts, set variable start_time2 to value given value
        nRowsForTrajectoryFile = T.x.shape[0]
        timeGrain = T.time[1]
        trial_length = (nRowsForTrajectoryFile - 1) * timeGrain
        
        start_time2 = 0
        stop_time = cue_time + trial_length 
        
        P.update_positions_psychopy(start_time2)
        pos_obj1 = P.objects[0].pos
        
        # keep track of which components have finished
        trialComponents = [P.objects[0], cue, noise]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        P.objects[0].finalx0 = -999 #this is just for recording the final position to the data file
        P.objects[0].finaly0 = -999
        P.objects[0].penultimatex0 = -999 #this is just for recording the penultimate position to the data file
        P.objects[0].penultimatey0 = -999 #this is just for recording the penultimate position to the data file
        P.objects[0].antepenultimatex0 = -999 #this is just for recording the antepenultimate position to the data file
        P.objects[0].antepenultimatey0 = -999 #this is just for recording the antepenultimate position to the data file
        P.objects[0].preantepenultimatex0 = -999 #this is just for recording the preantepenultimate position to the data file
        P.objects[0].preantepenultimatey0 = -999 #this is just for recording the preantepenultimate position to the data file
        time.sleep(.2)
        # reset timers
        t = 0
        nKeys = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        finishedCriticalStimuli = False
        
        # -------Run Routine "trial"-------
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=trialClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            timeSuspectGlobalClockAheadBy = tThisFlipGlobal - t
            globalClockAheadBy = trialClock.getLastResetTime() - logging.defaultClock.getLastResetTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            noise.draw( )
            o1.setAutoDraw(False)
            cue.setAutoDraw(False)
            # cueing of the targets phase at beginning of trial
            if t >= 0 and t <= cue_time:
                text.draw( )
                # simple iterator through frames
                # we are not using t, as it is not an integer (we could round though)
                iframe = iframe + 1     
                
                # update and draw objects
                P.update_positions_psychopy(start_time2)    
                for i2 in range(n_objects):
                    P.objects[i2].draw() 
                #move phase
            if t >= cue_time and t <= stop_time:
                event.clearEvents()
                text.draw( )
                for i2 in range(nTargets):
                    P.objects[i2].size = size_normal
                P.update_positions_psychopy(start_time2 + t-cue_time)
                for i2 in range(n_objects):
                    P.objects[i2].draw()
                numOfFinalFramesToRecord = 40
                tLastFramesStart = stop_time - numOfFinalFramesToRecord*(1.0/refreshRateObserved) #time after which should be the last few frames
                if t > tLastFramesStart: #make sure pick up the last frames
                    #Record the time of the very last frame that the stimuli are still on, by doing it over and over, so when the code ceases it will be correct.
                    P.objects[0].tLastFrame = -999 #dummy value which will hopefully get overwritten on win flip by the next line, over and over until the very last frame
                    win.timeOnFlip(P.objects[0], 'tLastFrame')  # set P.objects[0].tLastFrame to time at next scr refresh
                
                    #To record the penultimate position to the data file, on each pass-through, set it to final before setting final to current, so that
                    #at the last frame, final will be final
                    P.objects[0].preantepenultimatex0 = P.objects[0].antepenultimatex0
                    P.objects[0].preantepenultimatey0 = P.objects[0].antepenultimatey0
                    P.objects[0].antepenultimatex0 = P.objects[0].penultimatex0
                    P.objects[0].antepenultimatey0 = P.objects[0].penultimatey0               
                    P.objects[0].penultimatex0 = P.objects[0].finalx0 #You might think could rely on whatever the final x is, but this making sure recording the final x *drawn*
                    P.objects[0].penultimatey0 = P.objects[0].finaly0 #You might think could rely on whatever the final x is, but this making sure recording the final x *drawn*
                    #Record last x,y of object 0
                    P.objects[0].finalx0 = P.objects[0].pos[0] #You might think could rely on whatever the final x is, but this making sure recording the final x *drawn*
                    P.objects[0].finaly0 = P.objects[0].pos[1]
                    #But win.timeOnFlip won't work for recording a whole list or array of times because timeOnFlip can't append to an array.
                    #So for recording the interframe intervals of these times the conventional way is https://www.psychopy.org/general/timing/detectingFrameDrops.html
                    win.recordFrameIntervals = True
                    
                    #Another way is to record it directly myself
                    #https://github.com/psychopy/psychopy/blob/e45520f446697d5b5fad035fa0e4e50088ffa535/psychopy/visual/window.py
                    #That means calling my own function that gets the time of the flip.
                    #To get the time of the flip, it seems that the Psychopy functions use win._frameTime, which is set by logging.defaultClock.getTime() automatically
                    #Alternatively, I could probably use the automatic log of win._frameTimes, but then at the time of calling addData, to know 
                    #which are the critical frames, I'll need to record how long the array of win._frameTimes is at the last frame.
                    #win.callOnFlip(appendToLastFrameTimes, frameTime=win._frameTime
                    #win.timeOnFlip( lastFrameTimes[0], 
                
                
            # query phase
            if t > stop_time:
                text.draw( )
                keys = event.getKeys()
                for thisKey in keys:
                    if (thisKey == 'up'):
                        wasUp = 1
                        nKeys = nKeys + 1
                    elif (thisKey == 'down'): 
                        wasUp = 0
                        nKeys = nKeys + 1
                if not finishedCriticalStimuli:
                    lenFrameTimesUntilStimuliFinish = len(win._frameTimes)
                finishedCriticalStimuli = True
                win.recordFrameIntervals = False
                if 'dynamic' in thisStair.extraInfo:
                    if noise.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                        noise.frameNStart = frameN  # exact frame index
                        noise.tStart = t  # local t and not account for scr refresh
                        noise.tStartRefresh = tThisFlipGlobal  # on global time, but not as accurate as win.time so overrule it in next line with that
                        win.timeOnFlip(noise, 'tStartRefresh')  # time at next scr refresh
                        noise.status = STARTED #AOH
    
                    if noise.status == STARTED:
                        if noise._needBuild:
                            noise.buildNoise()
                        else:
                            noise.draw() #AOH
                            if (frameN-noise.frameNStart) % 2 == 0:
                                noise.updateNoise()
                else:
                    noise.draw( )
                if nKeys > 0:
                    thisStair.addResponse(wasUp)
                    expInfo['response'] = keys
                    expInfo['wasUp'] = wasUp
                    expInfo['condition'] = thisStair.extraInfo
                    if 'f' in track_filename:
                        expInfo['flipStatus'] = 'flipped'
                    else: expInfo['flipStatus'] = 'not_flipped'
                    expInfo['offset'] = thisIntensity
                    trialNumber = trialNumber + 1
                    expInfo['trialNumber'] = trialNumber
                    continueRoutine = False
                    for i2 in range(n_objects):
                        P.objects[i2].setAutoDraw(False)
                    thisExp.saveAsWideText(filename+'.csv', fileCollisionMethod='overwrite',delim='auto')
                    time.sleep(0.05) 
            
            # *o1* updates
            if P.objects[0].status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                P.objects[0].frameNStart = frameN  # exact frame index
                P.objects[0].tStart = t  # local t and not account for scr refresh
                win.timeOnFlip(P.objects[0], 'tStartRefresh')  # time at next scr refresh
               
            if cue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                cue.frameNStart = frameN  # exact frame index
                cue.tStart = t  # local t and not account for scr refresh
                cue.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(cue, 'tStartRefresh')  # time at next scr refresh
                cue.setAutoDraw(True)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # save mouse data
        win.blendmode = 'avg'
    
        if debug:
            print('globalClockAheadBy=',globalClockAheadBy)
            print('P.objects[0].tStartRefresh - globalClockAheadBy = ',P.objects[0].tStartRefresh - globalClockAheadBy)
        trials.addData('P.objects[0].started', P.objects[0].tStartRefresh - globalClockAheadBy)
        trials.addData('P.objects[0].tLastFrame', P.objects[0].tLastFrame - globalClockAheadBy)
        trials.addData('P.objects[0].penultimatex0',P.objects[0].penultimatex0)
        trials.addData('P.objects[0].penultimatey0',P.objects[0].penultimatey0)
        trials.addData('P.objects[0].antepenultimatex0',P.objects[0].antepenultimatex0)
        trials.addData('P.objects[0].antepenultimatey0',P.objects[0].antepenultimatey0)
        trials.addData('P.objects[0].preantepenultimatex0',P.objects[0].preantepenultimatex0)
        trials.addData('P.objects[0].preantepenultimatey0',P.objects[0].preantepenultimatey0)
        trials.addData('P.objects[0].finalx0',P.objects[0].finalx0)
        trials.addData('P.objects[0].finaly0',P.objects[0].finaly0)
        
    #Record final frametimes when critical stimuli on
        trials.addData('win.nDroppedFrames',win.nDroppedFrames)
        intervals_msec = pylab.array(win.frameIntervals) * 1000
        intervals_msec_last_thistrial = intervals_msec[-numOfFinalFramesToRecord:] #last numOfFinalFramesToRecord
        #intervals_msec = intervals_msec[numOfFinalFramesToRecord
        if debug:
            print('Frame times up to stimuli cessation according to win.frameIntervals=', intervals_msec_last_thistrial)
            print('In last',numOfFinalFramesToRecord,' frames across all trials ', win.nDroppedFrames, ' frames were dropped according to nDroppedFrames.')
        trials.addData('finalStimFrameTimes',intervals_msec_last_thistrial)
        worstLastFramesDur = max(intervals_msec_last_thistrial)
        trials.addData('worstLastFramesDur',worstLastFramesDur)
        #Count the number of timing hiccups and print out and save to log some information about them
        longFrameLimit = np.round(1000/refreshRateObserved*(1.0+frameTimeTolerance),2)
        idxsInterframeLong = np.where( np.array(intervals_msec_last_thistrial) > longFrameLimit ) [0] #frames that exceeded frameTimeTolerance of expected duration
        numCasesInterframeLong = len( idxsInterframeLong )
        trials.addData('timingHiccupsInLastFramesOfStimuli',numCasesInterframeLong)
        if numCasesInterframeLong >0 and (not demo):
           longFramesStr =  'ERROR,'+str(numCasesInterframeLong)+' frames were longer than '+str(longFrameLimit)+' ms'
           if demo: 
             longFramesStr += 'not printing them all because in demo mode'
           else:
               longFramesStr += ' apparently screen refreshes skipped, interframe durs were:'+\
                        str( np.around(  intervals_msec_last_thistrial[idxsInterframeLong] ,1  ) )+ ' and was these frames: '+ str(idxsInterframeLong)
           if longFramesStr != None:
                    logging.error( 'trialnum='+str(trials.thisN)+' '+longFramesStr )
                    if not demo:
                        flankingAlso=list()
                        for idx in idxsInterframeLong: #also print timing of one before and one after long frame
                            if idx-1>=0:
                                flankingAlso.append(idx-1)
                            else: flankingAlso.append(np.NaN)
                            flankingAlso.append(idx)
                            if idx+1<len(intervals_msec_last_thistrial):  flankingAlso.append(idx+1)
                            else: flankingAlso.append(np.NaN)
                        flankingAlso = np.array(flankingAlso)
                        flankingAlso = flankingAlso[np.logical_not(np.isnan(flankingAlso))]  #remove nan values
                        flankingAlso = flankingAlso.astype(np.integer) #cast as integers, so can use as subscripts
                        logging.info( 'flankers also='+str( np.around( intervals_msec_last_thistrial[flankingAlso], 1) )  ) #because this is not an essential error message, as previous one already indicates error
                        #As INFO,it won't fill up the console when console set to WARNING or higher
        # store data for trials (TrialHandler)
        trials.addData('cue.started', cue.tStartRefresh)
        trials.addData('cue.stopped', cue.tStopRefresh)
        if debug and (practice_trials.thisN==2): #Quit and Plot frame intervals
            win.close()
            plotFrameIntervals(intervals_msec)
            core.quit()
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
    
# completed 1 repeats of 'trials'


# ------Prepare to start Routine "bye"-------
continueRoutine = True
routineTimer.add(2.000000)
# update component parameters for each repeat
# keep track of which components have finished
byeComponents = [text_4]
for thisComponent in byeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
byeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "bye"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = byeClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=byeClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_4* updates
    if text_4.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_4.frameNStart = frameN  # exact frame index
        text_4.tStart = t  # local t and not account for scr refresh
        text_4.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_4, 'tStartRefresh')  # time at next scr refresh
        text_4.setAutoDraw(True)
    if text_4.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > text_4.tStartRefresh + 2.0-frameTolerance:
            # keep track of stop time/frame for later
            text_4.tStop = t  # not accounting for scr refresh
            text_4.frameNStop = frameN  # exact frame index
            win.timeOnFlip(text_4, 'tStopRefresh')  # time at next scr refresh
            text_4.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in byeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "bye"-------
for thisComponent in byeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
