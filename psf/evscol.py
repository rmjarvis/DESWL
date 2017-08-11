import matplotlib.pyplot as plt
import numpy

# Make fake data.  The real thing will use the real ngmix011 data.
cols = numpy.arange(20,2030,5)
e1 = numpy.random.uniform( -0.001, 0.001, len(cols) )
e2 = numpy.random.uniform( -0.001, 0.001, len(cols) )

# Mock up the slop at each end
e1[:5] -= (45-cols[:5])*0.006/25.
e1[-5:] -= (cols[-5:]-2005)*0.006/25.

# Make 3 axes, where the first one spans the whole row.
# cf. http://matplotlib.org/users/gridspec.html
plt.style.use('supermongo')
fig = plt.figure()
ax1 = plt.subplot2grid( (2,2), (0,0), colspan=2 )
ax2 = plt.subplot2grid( (2,2), (1,0) )
ax3 = plt.subplot2grid( (2,2), (1,1) )

# Plot data on first three axes:
ax1.scatter(cols, e1, marker='o', color='blue', s=1.0, label=r'$\langle e_1 \rangle$')
ax1.scatter(cols, e2, marker='o', color='green', s=1.0, label=r'$\langle e_2 \rangle$')

ax2.scatter(cols, e1, marker='o', color='blue', s=1.5)
ax2.scatter(cols, e2, marker='o', color='green', s=1.5)

ax3.scatter(cols, e1, marker='o', color='blue', s=1.5)
ax3.scatter(cols, e2, marker='o', color='green', s=1.5)

# Draw the legend only once
ax1.legend(loc=(0.75, 0.08), fontsize=15)

# Limit the axis ranges
ax1.set_xlim(0, 2047)    # Show the whole range, but no points in first or last 20.
ax2.set_xlim(0, 110)     # Note: go 10 past where the last point will be drawn.
ax3.set_xlim(1937,2047)  # Likewise, start 10 before first point.

ymin = -0.006
ymax = 0.003
ax1.set_ylim(ymin, ymax)
ax2.set_ylim(ymin, ymax)
ax3.set_ylim(ymin, ymax)

# Make ax2, ax3 look like a single plot with a broken x-axis
ax2.spines['right'].set_visible(False) # Hide the right spine
ax2.yaxis.tick_left()                  # Only put ticks on the left.

ax3.spines['left'].set_visible(False)  # Hide the left spine
ax3.yaxis.tick_right()                 # Only put ticks on the right.
ax3.yaxis.set_ticklabels([])           # Don't label the y ticks.

# It wants to label ever 0.001, but I think it looks nicer every 0.002.
ax1.yaxis.set_ticks( numpy.arange(-0.006, 0.004, 0.002) )
ax2.yaxis.set_ticks( numpy.arange(-0.006, 0.004, 0.002) )
ax3.yaxis.set_ticks( numpy.arange(-0.006, 0.004, 0.002) )

# Make little diagonal cuts to make the discontinuity clearer:
# cf. http://stackoverflow.com/questions/5656798/python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis
d = .03 # how big to make the diagonal lines in axes coordinates
m = 0.7  # slope of the cut lines
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
ax2.plot((1-m*d,1+m*d),(-d,+d), **kwargs) # top-left diagonal
ax2.plot((1-m*d,1+m*d),(1-d,1+d), **kwargs) # bottom-left diagonal

kwargs.update(transform=ax3.transAxes) # switch to the bottom axes
ax3.plot((-m*d,m*d),(-d,+d), **kwargs) # top-right diagonal
ax3.plot((-m*d,m*d),(1-d,1+d), **kwargs) # bottom-right diagonal

# Squeeze them a bit closer together
plt.subplots_adjust(wspace=0.06)

# Put the y label on ax1, ax2
ax1.set_ylabel(r'$\langle e \rangle$')
ax2.set_ylabel(r'$\langle e \rangle$')
# Make a little more room for the label.
plt.subplots_adjust(left=0.16)

# For the x axis, we want it to be centered.  Easiest to just do this by hand.
fig.text(0.5, 0.04, 'X Position on Chip', ha='center', va='center')
# Make room.
plt.subplots_adjust(bottom=0.12)

# Shade in a grey region where the chips are masked
ax1.fill_between([0,15],[ymin,ymin], [ymax,ymax], color='LightGrey')
ax1.fill_between([2032,2047],[ymin,ymin], [ymax,ymax], color='LightGrey')
ax2.fill_between([0,15],[ymin,ymin], [ymax,ymax], color='LightGrey')
ax3.fill_between([2032,2047],[ymin,ymin], [ymax,ymax], color='LightGrey')

plt.savefig('evscol.eps')  # This is usually the best for putting in the paper.
plt.savefig('evscol.png')  # This is usually quicker to view to check how things look.
