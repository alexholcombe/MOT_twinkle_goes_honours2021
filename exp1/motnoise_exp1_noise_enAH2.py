#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2021.1.2),
    on Tue Mar 23 14:13:24 2021
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding
import time
import copy
from psychopy.hardware import keyboard


# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2021.1.2'
expName = 'noiseMot_exp1_noise'  # from the Builder filename that created this script
expInfo = {'motbox_path': '../', 'participant': '999', 'protocol': '1'}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='/Users/alexh/Documents/Supervision_rec_ltrs_recommendations_referee/SupervisionAndRecLtrs/honours_Dissertatns/2021honoursSupervisionAlex/motNoise_forAlex/exp1/motnoise_exp1_noise_en.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.DATA)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# Setup the Window
win = visual.Window(
    size=[1920, 1200], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[0.004,0.004,0.004], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='norm')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "welcome"
welcomeClock = core.Clock()
text_2 = visual.TextStim(win=win, name='text_2',
    text='In this experiment, your task is to track moving objects. In each trial, there will be 8 static objects, while 4 of them will blink. After 2 seconds, all objects became indistinguishible and began to move. Your task is to track the originally highlighted objects.\nAfter the movement ends, your task is to click on the originally highlighted objects\n\nPress space to continue..',
    font='Arial',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
key_resp_3 = keyboard.Keyboard()

# Initialize components for Routine "intro_block"
intro_blockClock = core.Clock()
text_3 = visual.TextStim(win=win, name='text_3',
    text='First, there will be training part (12 trials). You will see three types of difficulty (from easiest to most difficult)\n\nPress space',
    font='Arial',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
#gabor_preview_pth = 'gabors/gabor_contr{:.1f}.png'.format(float(expInfo['t_contr']))

prot_file = 'protocols/P%03d_noise.csv' % int(expInfo['protocol'])
key_resp_4 = keyboard.Keyboard()

# Initialize components for Routine "fix"
fixClock = core.Clock()
text = visual.TextStim(win=win, name='text',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
image = visual.ImageStim(
    win=win,
    name='image', 
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)

# Initialize components for Routine "trial"
trialClock = core.Clock()

sys.path.insert(0, expInfo['motbox_path']) # to load module from neighbouring folder
import motbox

nTargets = 4 # number of targets (that will be highlighted in cue phase)
n_objects = 8 # number of objects (could be loaded automatically from file

# configure experiment
trial_length = 6 # duration of trial in seconds
cue_time = 2     # duration of cueing phase in seconds

# hide mouse cursor
win.mouseVisible = False

mouse = event.Mouse(win=win)
x, y = [None, None]
mouse.mouseClock = core.Clock()
p1 = visual.Polygon(
    win=win, name='p1',units='pix', 
    edges=50, size=(64, 64),
    ori=0, pos=(10000, 10000),
    lineWidth=1,     colorSpace='rgb',  lineColor=[1,1,1], fillColor=None,
    opacity=1, depth=-3.0, interpolate=True)

# Initialize components for Routine "before_trials"
before_trialsClock = core.Clock()
text_5 = visual.TextStim(win=win, name='text_5',
    text='Now the actual experiment begins\n\n\n\n\nPress space',
    font='Arial',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
key_resp_2 = keyboard.Keyboard()

# Initialize components for Routine "fix"
fixClock = core.Clock()
text = visual.TextStim(win=win, name='text',
    text='+',
    font='Arial',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);
image = visual.ImageStim(
    win=win,
    name='image', 
    image='sin', mask=None,
    ori=0, pos=(0, 0), size=(2, 2),
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)

# Initialize components for Routine "trial"
trialClock = core.Clock()


sys.path.insert(0, expInfo['motbox_path']) # to load module from neighbouring folder
import motbox
import time
import copy


nTargets = 4 # number of targets (that will be highlighted in cue phase)
n_objects = 8 # number of objects (could be loaded automatically from file

# configure experiment
trial_length = 6 # duration of trial in seconds
cue_time = 2     # duration of cueing phase in seconds

# hide mouse cursor
win.mouseVisible = False

mouse = event.Mouse(win=win)
x, y = [None, None]
mouse.mouseClock = core.Clock()
p1 = visual.Polygon(
    win=win, name='p1',units='pix', 
    edges=50, size=(64, 64),
    ori=0, pos=(10000, 10000),
    lineWidth=1,     colorSpace='rgb',  lineColor=[1,1,1], fillColor=None,
    opacity=1, depth=-3.0, interpolate=True)

# Initialize components for Routine "bye"
byeClock = core.Clock()
text_4 = visual.TextStim(win=win, name='text_4',
    text='Thanks!',
    font='Arial',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

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
        win.timeOnFlip(text_2, 'tStartRefresh')  # time at next scr refresh
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
thisExp.addData('text_2.started', text_2.tStartRefresh)
thisExp.addData('text_2.stopped', text_2.tStopRefresh)
# the Routine "welcome" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

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
thisExp.addData('text_3.started', text_3.tStartRefresh)
thisExp.addData('text_3.stopped', text_3.tStopRefresh)
win.blendMode = 'avg'
# the Routine "intro_block" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
practice_trials = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('protocols' + os.path.sep + 'P_practice_noise.csv'),
    seed=None, name='practice_trials')
thisExp.addLoop(practice_trials)  # add the loop to the experiment
thisPractice_trial = practice_trials.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisPractice_trial.rgb)
if thisPractice_trial != None:
    for paramName in thisPractice_trial:
        exec('{} = thisPractice_trial[paramName]'.format(paramName))

for thisPractice_trial in practice_trials:
    currentLoop = practice_trials
    # abbreviate parameter names if possible (e.g. rgb = thisPractice_trial.rgb)
    if thisPractice_trial != None:
        for paramName in thisPractice_trial:
            exec('{} = thisPractice_trial[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "fix"-------
    continueRoutine = True
    routineTimer.add(0.500000)
    # update component parameters for each repeat
    noise_backg_pth = 'noise/noise%03d.png' % noise_id
    win.blendMode = 'add'
    image.setImage(noise_backg_pth)
    # keep track of which components have finished
    fixComponents = [text, image]
    for thisComponent in fixComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    fixClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "fix"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = fixClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=fixClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text* updates
        if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text.frameNStart = frameN  # exact frame index
            text.tStart = t  # local t and not account for scr refresh
            text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
            text.setAutoDraw(True)
        if text.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > text.tStartRefresh + 0.5-frameTolerance:
                # keep track of stop time/frame for later
                text.tStop = t  # not accounting for scr refresh
                text.frameNStop = frameN  # exact frame index
                win.timeOnFlip(text, 'tStopRefresh')  # time at next scr refresh
                text.setAutoDraw(False)
        win.blendMode = 'avg'
        
        # *image* updates
        if image.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            image.frameNStart = frameN  # exact frame index
            image.tStart = t  # local t and not account for scr refresh
            image.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(image, 'tStartRefresh')  # time at next scr refresh
            image.setAutoDraw(True)
        if image.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > image.tStartRefresh + 0.5-frameTolerance:
                # keep track of stop time/frame for later
                image.tStop = t  # not accounting for scr refresh
                image.frameNStop = frameN  # exact frame index
                win.timeOnFlip(image, 'tStopRefresh')  # time at next scr refresh
                image.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in fixComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "fix"-------
    for thisComponent in fixComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    practice_trials.addData('text.started', text.tStartRefresh)
    practice_trials.addData('text.stopped', text.tStopRefresh)
    win.blendMode = 'avg'
    practice_trials.addData('image.started', image.tStartRefresh)
    practice_trials.addData('image.stopped', image.tStopRefresh)
    
    # ------Prepare to start Routine "trial"-------
    continueRoutine = True
    # update component parameters for each repeat
    # probably redudndant
    win.mouseVisible = False
    
    # load gabor based on the used contrast
    gabor_name_pth = 'gabors/gabor_contr_id_%01d.png' % t_contr
    
    # noise needs to be created each routine, as I run into some problems when reusing old objects
    noise_backg = visual.ImageStim(
        win=win, name='image',
        image=noise_backg_pth, mask=None,
        ori=0, pos=(0, 0), size=(2, 2),
        color=[1,1,1], colorSpace='rgb', opacity=1,
        flipHoriz=False, flipVert=False,
        texRes=128, interpolate=True, depth=-3.0)
    
    # same for the objects. This is a placeholder for one object
    o1 = visual.Circle(
        win=win, name='o1',units='pix', 
        pos=(10000, 10000), radius=64,
        lineColor=[1,1,1], fillColor=[1,1,1], color=[1,1,1],
        colorSpace='rgb', opacity=1)
#    o1 = visual.ImageStim(
#        win=win, name='o1',units='pix', 
#        image=gabor_name_pth, mask=None,
#        ori=0, pos=(10000, 10000), size=(64, 64),
#        color=[1,1,1], colorSpace='rgb', opacity=1,
#        flipHoriz=False, flipVert=False,
#        texRes=128, interpolate=True, depth=-1.0)
    
    # and circular highlighting
    # this is related to current experiment (see poster)
    p1 = visual.Polygon(
        win=win, name='p1',units='pix', 
        edges=50, size=(64, 64),
        ori=0, pos=(10000, 10000),
        lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
        fillColor=None, fillColorSpace='rgb',
        opacity=1, depth=-3.0, interpolate=True)

    # load trajectory from file (this was generated using motrack package
    track_filename = 'trajectories/T%03d.csv' % trajectory_id
    
    # some intitial setup
    clicked_name_correct = [];
    mouse_x_correct = [];
    mouse_y_correct = [];
    mouse_time_correct = [];
    
    # initialize motbox trajectory
    T = motbox.Track()
    T.load_from_csv(track_filename, delim = ',') # load it from file
    
    # scale it with given pixels per inch. This needs to be set for proper values
    T.scale(30) # px per inch 
    
    # Puppeteer handles the actual moving of objects
    P = motbox.Puppeteer()
    P.track = T
    # this creates multiple copies from the template
    P.clone_template_psychopy(o1, n_objects)
    
    # create polygons highlighting targets
    for x_obj3 in P.objects:
        x_obj3.pressed   = False
        x_obj3.mmark      = copy.copy(p1)
        x_obj3.mmark.name = "p1_{}".format(x_obj3.name)
    
    # this is related to current experiment, in half of trials, objects were highlighted in query phase
    shouldMark = mark_type == 'mark'
    nClicks = 0
    
    win.blendMode = 'avg' #'add' #set proper blending mode AOH  Average (avg) places the new stimulus over the old one with a transparency given by its opacity.
    iframe = 1 # this is iterator through frames
    
    size_blink  = 8       # how many frames are there between size change 
    size_big_coef    = 2  # how large are objects when their size is increased
    size_normal_coef = 1  # normal size coeficient, should be 1
    size_normal = P.objects[0].size
    size_coef = size_normal_coef 
    
    # set start time
    # in this particular study, we are not randomizing starts of the trials
    # to randomize starts, set variable start_time2 to value given value
    start_time2 = 0
    stop_time = cue_time + trial_length 
    
    # mark
    first_mark = 1
    #o1.setImage(gabor_name_pth) #AOH
    
    # setup some python lists for storing info about the mouse
    mouse.x = []
    mouse.y = []
    mouse.leftButton = []
    mouse.midButton = []
    mouse.rightButton = []
    mouse.time = []
    mouse.clicked_name = []
    gotValidClick = False  # until a click is received
    mouse.mouseClock.reset()
    # keep track of which components have finished
    trialComponents = [o1, mouse, p1]
    for thisComponent in trialComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "trial"-------
    while continueRoutine:
        # get current time
        t = trialClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=trialClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        noise_backg.setAutoDraw(False) #AOH
        o1.setAutoDraw(False)
        p1.setAutoDraw(False)
        
        # cue phase
        if t >= 0 and t <= cue_time:
            
            # simple iterator through frames
            # we are not using t, as it is not an integer (we could round though)
            iframe = iframe + 1     
            
            # whether we should change size
            # there should be color changing or something else, but it was difficult to 
            # get proper color schemes with blending
            if (iframe % size_blink) == 0:
                if size_coef == size_normal_coef:
                    size_coef = size_big_coef
                else:
                    size_coef = size_normal_coef
            # change size for targets
            for i1 in range(nTargets):
                P.objects[i1].size = size_normal * size_coef
            
            # update and draw objects
            P.update_positions_psychopy(start_time2)    
            for i2 in range(n_objects):
                P.objects[i2].draw() 
            
        # move phase
        if t >= cue_time and t <= stop_time:
            # normalize sizes 
            for i2 in range(nTargets):
                P.objects[i2].size = size_normal
            # update and draw objects
            P.update_positions_psychopy(start_time2 + t-cue_time)
            for i2 in range(n_objects):
                P.objects[i2].draw() 
            
            
        # query phase
        if t > stop_time:
            # this is used only in trials with target highlighting
            # basically it draws a circle around the objects
            if (first_mark and shouldMark):
                first_mark = 0;
                for x_obj2 in P.objects:
                    x_obj2.mmark.pos = x_obj2.pos
                    x_obj2.mmark.setAutoDraw(True)
            # repeat, until participant clicked on 4 objects
            if (nClicks < nTargets):
                win.mouseVisible = True
                #nClicks = nClicks + 1 #AOH to use this once understand query code better
                # this code handles that we can't click on same object multiple times
                for x_obj2 in P.objects:
                    if mouse.isPressedIn(x_obj2, buttons=[0]) and x_obj2.pressed == False:
                        x, y = mouse.getPos()
                        x_obj2.pressed  = True
                        nClicks = nClicks + 1 #original code
                        x_obj2.size = x_obj2.size * 2 # clicked object increases its size
                        
                        # add information about the clicked objects
                        clicked_name_correct.append(x_obj2.name)
                        mouse_x_correct.append(x)
                        mouse_y_correct.append(y)
                        mouse_time_correct.append(trialClock.getTime())
                # update and draw objects (note that 'update' keeps object at the same locations
                P.update_positions_psychopy(start_time2 + stop_time-cue_time)
                for i2 in range(n_objects):
                    P.objects[i2].draw() 
            else: # all four targets were clicked, ending routine
                win.mouseVisible = False
                continueRoutine = False
                for i2 in range(n_objects):
                    P.objects[i2].mmark.setAutoDraw(False)
                    P.objects[i2].setAutoDraw(False)
                noise_backg.setAutoDraw(False)
                time.sleep(0.3) 
        
        # *o1* updates
        if o1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            o1.frameNStart = frameN  # exact frame index
            o1.tStart = t  # local t and not account for scr refresh
            o1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(o1, 'tStartRefresh')  # time at next scr refresh
            o1.setAutoDraw(True)
        # *mouse* updates
        if mouse.status == NOT_STARTED and t >= 4-frameTolerance:
            # keep track of start time/frame for later
            mouse.frameNStart = frameN  # exact frame index
            mouse.tStart = t  # local t and not account for scr refresh
            mouse.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(mouse, 'tStartRefresh')  # time at next scr refresh
            mouse.status = STARTED
            prevButtonState = mouse.getPressed()  # if button is down already this ISN'T a new click
        if mouse.status == STARTED:  # only update if started and not finished!
            buttons = mouse.getPressed()
            if buttons != prevButtonState:  # button state changed?
                prevButtonState = buttons
                if sum(buttons) > 0:  # state changed to a new click
                    #gotValidClick = True #AOH consider any click, regardless of position, to be valid
                    # check if the mouse was inside our 'clickable' objects
                    for obj in [P.objects[0],P.objects[1],P.objects[2],P.objects[3],P.objects[4],P.objects[5],P.objects[6],P.objects[7]]:
                        if obj.contains(mouse):
                            gotValidClick = True
                            gotClickOnObject = True #AOH to rely on this, once he understands the two sets of code
                            mouse.clicked_name.append(obj.name)
                    x, y = mouse.getPos()
                    mouse.x.append(x)
                    mouse.y.append(y)
                    buttons = mouse.getPressed()
                    mouse.leftButton.append(buttons[0])
                    mouse.midButton.append(buttons[1])
                    mouse.rightButton.append(buttons[2])
                    mouse.time.append(mouse.mouseClock.getTime())
        
        # *p1* updates
        if p1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            p1.frameNStart = frameN  # exact frame index
            p1.tStart = t  # local t and not account for scr refresh
            p1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(p1, 'tStartRefresh')  # time at next scr refresh
            p1.setAutoDraw(True)
        
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
    mouse.clicked_name = clicked_name_correct
    mouse.x = mouse_x_correct
    mouse.y = mouse_y_correct
    mouse.time = mouse_time_correct
    
    practice_trials.addData('o1.started', o1.tStartRefresh)
    practice_trials.addData('o1.stopped', o1.tStopRefresh)
    # store data for practice_trials (TrialHandler)
    practice_trials.addData('mouse.x', mouse.x)
    practice_trials.addData('mouse.y', mouse.y)
    practice_trials.addData('mouse.leftButton', mouse.leftButton)
    practice_trials.addData('mouse.midButton', mouse.midButton)
    practice_trials.addData('mouse.rightButton', mouse.rightButton)
    practice_trials.addData('mouse.time', mouse.time)
    practice_trials.addData('mouse.clicked_name', mouse.clicked_name)
    practice_trials.addData('mouse.started', mouse.tStart)
    practice_trials.addData('mouse.stopped', mouse.tStop)
    practice_trials.addData('p1.started', p1.tStartRefresh)
    practice_trials.addData('p1.stopped', p1.tStopRefresh)
    # the Routine "trial" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# completed 1 repeats of 'practice_trials'


# ------Prepare to start Routine "before_trials"-------
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
        theseKeys = key_resp_2.getKeys(keyList=['space'], waitRelease=False)
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
thisExp.addData('text_5.started', text_5.tStartRefresh)
thisExp.addData('text_5.stopped', text_5.tStopRefresh)
# the Routine "before_trials" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
trials = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions(prot_file),
    seed=None, name='trials')
thisExp.addLoop(trials)  # add the loop to the experiment
thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
if thisTrial != None:
    for paramName in thisTrial:
        exec('{} = thisTrial[paramName]'.format(paramName))

for thisTrial in trials:
    currentLoop = trials
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            exec('{} = thisTrial[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "fix"-------
    continueRoutine = True
    routineTimer.add(0.500000)
    # update component parameters for each repeat
    noise_backg_pth = 'noise/noise%03d.png' % noise_id
    win.blendMode = 'add'
    image.setImage(noise_backg_pth)
    # keep track of which components have finished
    fixComponents = [text, image]
    for thisComponent in fixComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    fixClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "fix"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = fixClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=fixClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text* updates
        if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text.frameNStart = frameN  # exact frame index
            text.tStart = t  # local t and not account for scr refresh
            text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
            text.setAutoDraw(True)
        if text.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > text.tStartRefresh + 0.5-frameTolerance:
                # keep track of stop time/frame for later
                text.tStop = t  # not accounting for scr refresh
                text.frameNStop = frameN  # exact frame index
                win.timeOnFlip(text, 'tStopRefresh')  # time at next scr refresh
                text.setAutoDraw(False)
        win.blendMode = 'avg'
        
        # *image* updates
        if image.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            image.frameNStart = frameN  # exact frame index
            image.tStart = t  # local t and not account for scr refresh
            image.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(image, 'tStartRefresh')  # time at next scr refresh
            image.setAutoDraw(True)
        if image.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > image.tStartRefresh + 0.5-frameTolerance:
                # keep track of stop time/frame for later
                image.tStop = t  # not accounting for scr refresh
                image.frameNStop = frameN  # exact frame index
                win.timeOnFlip(image, 'tStopRefresh')  # time at next scr refresh
                image.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in fixComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "fix"-------
    for thisComponent in fixComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    trials.addData('text.started', text.tStartRefresh)
    trials.addData('text.stopped', text.tStopRefresh)
    win.blendMode = 'avg'
    trials.addData('image.started', image.tStartRefresh)
    trials.addData('image.stopped', image.tStopRefresh)
    
    # ------Prepare to start Routine "trial"-------
    continueRoutine = True
    # update component parameters for each repeat
    # probably redudndant
    win.mouseVisible = False
    
    # load gabor based on the used contrast
    gabor_name_pth = 'gabors/gabor_contr_id_%01d.png' % t_contr
    
    # noise needs to be created each routine, I run into some problems, when reusing old objects
    noise_backg = visual.ImageStim(
        win=win, name='image',
        image=noise_backg_pth, mask=None,
        ori=0, pos=(0, 0), size=(2, 2),
        color=[1,1,1], colorSpace='rgb', opacity=1,
        flipHoriz=False, flipVert=False,
        texRes=128, interpolate=True, depth=-3.0)
    
    # same for the objects. This is a placeholder for one object
    o1 = visual.ImageStim(
        win=win, name='o1',units='pix', 
        image=gabor_name_pth, mask=None,
        ori=0, pos=(10000, 10000), size=(64, 64),
        color=[1,1,1], colorSpace='rgb', opacity=1,
        flipHoriz=False, flipVert=False,
        texRes=128, interpolate=True, depth=-1.0)
    
    # and circular highlighting
    # this is related to current experiment (see poster)
    p1 = visual.Polygon(
        win=win, name='p1',units='pix', 
        edges=50, size=(64, 64),
        ori=0, pos=(10000, 10000),
        lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
        fillColor=None, fillColorSpace='rgb',
        opacity=1, depth=-3.0, interpolate=True)
    
    # load trajectory from file (this was generated using the motrack R package)
    track_filename = 'trajectories/T%03d.csv' % trajectory_id
    
    # some intitial setup
    clicked_name_correct = [];
    mouse_x_correct = [];
    mouse_y_correct = [];
    mouse_time_correct = [];
    
    # initialize motbox trajectory
    T = motbox.Track()
    T.load_from_csv(track_filename, delim = ',') # load it from file
    
    # scale it with given pixels per inch. This needs to be set for proper values
    T.scale(30) # px per inch 
    
    # Puppeteer handles the actual moving of objects
    P = motbox.Puppeteer()
    P.track = T
    # this creates multiple copies from the template
    P.clone_template_psychopy(o1, n_objects)
    
    # create polygons highlighting targets
    for x_obj3 in P.objects:
        x_obj3.pressed   = False
        x_obj3.mmark      = copy.copy(p1)
        x_obj3.mmark.name = "p1_{}".format(x_obj3.name)
    
    # this is related to current experiment, in half of trials, objects were highlighted in query phase
    shouldMark = mark_type == 'mark'
    nClicks = 0
    
    win.blendMode = 'add' # set proper blending mode
    iframe = 1 # this is iterator through frames
    
    size_blink  = 8       # how many frames are there between size change 
    size_big_coef    = 2  # how large are objects when their size is increased
    size_normal_coef = 1  # normal size coeficient, should be 1
    size_normal = P.objects[0].size
    size_coef = size_normal_coef 
    
    # set start time
    # in this particular study, we are not randomizing starts of the trials
    # to randomize starts, set variable start_time2 to value given value
    
    start_time2 = 0
    stop_time = cue_time + trial_length 
    
    # mark
    first_mark = 1
    o1.setImage(gabor_name_pth)
    # setup some python lists for storing info about the mouse
    mouse.x = []
    mouse.y = []
    mouse.leftButton = []
    mouse.midButton = []
    mouse.rightButton = []
    mouse.time = []
    mouse.clicked_name = []
    gotValidClick = False  # until a click is received
    mouse.mouseClock.reset()
    # keep track of which components have finished
    trialComponents = [o1, mouse, p1]
    for thisComponent in trialComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "trial"-------
    while continueRoutine:
        # get current time
        t = trialClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=trialClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        noise_backg.setAutoDraw(True)
        o1.setAutoDraw(False)
        p1.setAutoDraw(False)
        
        # cue phase
        if t >= 0 and t <= cue_time:
            
            # simple iterator through frames
            # we are not using t, as it is not an integer (we could round though)
            iframe = iframe + 1     
            
            # whether we should change size
            # there should be color changing or something else, but it was difficult to 
            # get proper color schemes with blending
            if (iframe % size_blink) == 0:
                if size_coef == size_normal_coef:
                    size_coef = size_big_coef
                else:
                    size_coef = size_normal_coef
            # change size for targets
            for i1 in range(nTargets):
                P.objects[i1].size = size_normal * size_coef
            
            # update and draw objects
            P.update_positions_psychopy(start_time2)    
            for i2 in range(n_objects):
                P.objects[i2].draw() 
            
        # move phase
        if t >= cue_time and t <= stop_time:
            # normalize sizes 
            for i2 in range(nTargets):
                P.objects[i2].size = size_normal
            # update and draw objects
            P.update_positions_psychopy(start_time2 + t-cue_time)
            for i2 in range(n_objects):
                P.objects[i2].draw() 
            
        # query phase
        if t > stop_time:
            # this is used only in trials with target highlighting
            # basically, it draws a circle around the objects
            if (first_mark and shouldMark):
                first_mark = 0;
                for x_obj2 in P.objects:
                    x_obj2.mmark.pos = x_obj2.pos
                    x_obj2.mmark.setAutoDraw(True)
            # repeat, until participant clicked on 4 objects
            if (nClicks < nTargets):
                # we need mouse for that
                win.mouseVisible = True
                # this code handles that we can't click on same object multiple times
                for x_obj2 in P.objects:
                    if mouse.isPressedIn(x_obj2, buttons=[0]) and x_obj2.pressed == False:
                        x, y = mouse.getPos()
        
                        x_obj2.pressed  = True
                        nClicks = nClicks + 1
                        x_obj2.size = x_obj2.size * 2 # clicked object increases its size
                        
                        # add information about the clicked objects
                        clicked_name_correct.append(x_obj2.name)
                        mouse_x_correct.append(x)
                        mouse_y_correct.append(y)
                        mouse_time_correct.append(trialClock.getTime())
                # update and draw objects (note that 'update' keeps object at the same locations
                P.update_positions_psychopy(start_time2 + stop_time-cue_time)
                for i2 in range(n_objects):
                    P.objects[i2].draw() 
            
            else: # all four targets were clicked, ending routine
                win.mouseVisible = False
                continueRoutine = False
                for i2 in range(n_objects):
                    P.objects[i2].mmark.setAutoDraw(False)
                    P.objects[i2].setAutoDraw(False)
                noise_backg.setAutoDraw(False)
                time.sleep(0.3) 
        
        # *o1* updates
        if o1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            o1.frameNStart = frameN  # exact frame index
            o1.tStart = t  # local t and not account for scr refresh
            o1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(o1, 'tStartRefresh')  # time at next scr refresh
            o1.setAutoDraw(True)
        # *mouse* updates
        if mouse.status == NOT_STARTED and t >= 4-frameTolerance:
            # keep track of start time/frame for later
            mouse.frameNStart = frameN  # exact frame index
            mouse.tStart = t  # local t and not account for scr refresh
            mouse.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(mouse, 'tStartRefresh')  # time at next scr refresh
            mouse.status = STARTED
            prevButtonState = mouse.getPressed()  # if button is down already this ISN'T a new click
        if mouse.status == STARTED:  # only update if started and not finished!
            buttons = mouse.getPressed()
            if buttons != prevButtonState:  # button state changed?
                prevButtonState = buttons
                if sum(buttons) > 0:  # state changed to a new click
                    gotValidClick = True #AOH any click is valid, because task is to report target position
                    # check if the mouse was inside our 'clickable' objects
                    for obj in [P.objects[0],P.objects[1],P.objects[2],P.objects[3],P.objects[4],P.objects[5],P.objects[6],P.objects[7]]:
                        if obj.contains(mouse):
                            gotClickOnObject = True
                            mouse.clicked_name.append(obj.name)
                    x, y = mouse.getPos()
                    mouse.x.append(x)
                    mouse.y.append(y)
                    buttons = mouse.getPressed()
                    mouse.leftButton.append(buttons[0])
                    mouse.midButton.append(buttons[1])
                    mouse.rightButton.append(buttons[2])
                    mouse.time.append(mouse.mouseClock.getTime())
        
        # *p1* updates
        if p1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            p1.frameNStart = frameN  # exact frame index
            p1.tStart = t  # local t and not account for scr refresh
            p1.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(p1, 'tStartRefresh')  # time at next scr refresh
            p1.setAutoDraw(True)
        
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
    mouse.clicked_name = clicked_name_correct
    mouse.x = mouse_x_correct
    mouse.y = mouse_y_correct
    mouse.time = mouse_time_correct
    
    trials.addData('o1.started', o1.tStartRefresh)
    trials.addData('o1.stopped', o1.tStopRefresh)
    # store data for trials (TrialHandler)
    trials.addData('mouse.x', mouse.x)
    trials.addData('mouse.y', mouse.y)
    trials.addData('mouse.leftButton', mouse.leftButton)
    trials.addData('mouse.midButton', mouse.midButton)
    trials.addData('mouse.rightButton', mouse.rightButton)
    trials.addData('mouse.time', mouse.time)
    trials.addData('mouse.clicked_name', mouse.clicked_name)
    trials.addData('mouse.started', mouse.tStart)
    trials.addData('mouse.stopped', mouse.tStop)
    trials.addData('p1.started', p1.tStartRefresh)
    trials.addData('p1.stopped', p1.tStopRefresh)
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
thisExp.addData('text_4.started', text_4.tStartRefresh)
thisExp.addData('text_4.stopped', text_4.tStopRefresh)

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
