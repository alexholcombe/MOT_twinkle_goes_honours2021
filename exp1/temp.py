noise = visual.NoiseStim(
    win=win, name='noise',
    mask=None,
    ori=0.0, pos=(0, 0), size=(2048,2048),
    sf=None,
    color=[1,1,1], colorSpace='rgb', opacity=None, blendmode='avg', contrast=1.0,
    texRes=128, filter=None,
    noiseType='Uniform', noiseElementSize= (8, 8),
    interpolate=False, depth=0.0, units = 'pix')