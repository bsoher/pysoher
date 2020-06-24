"""This module contains the PlotPanel class that encapsulates a matplotlib
canvas within a wxpython widget.

The PlotPanel has a Figure and a Canvas and 'n' Axes. The user defines
the number of axes on Init and this number cannot be changed thereafter.
However, the user can select which of the n axes are actually displayed
in the Figure.

Axes are specified on Init because the zoom and reference cursors
need an axes to attach to to init properly.

on_size events simply set a flag, and the actual resizing of the figure is
triggered by an Idle event.

The mouse can be used to zoom in on the plot either across a 'span' or within
a 'box'. Mouse motion and button clicks can be bound to user defined methods.

PlotPanel Functionality
--------------------------------------------------
left mouse - If zoom mode is 'span', click and drag zooms the figure.
             A span is selected along the x-axis. On release, the
             axes xlim is adjusted accordingly. If zoom mode is 'box', then
             a zoom box is drawn during click and drag and figure is
             zoomed on both x- and y-axes upon release.

left mouse - click in place, un-zooms the figure to maximum x-data or
             x-data and y-data bounds.

right mouse - If reference mode is True/On, then click and drag will draw
              a span selector in the canvas that persists after release.

middle mouse - (or scroll button click), if do_middle_select_event and/or
               do_middle_motion_event are True then select, release and
               motion events are returned for these uses of the middle
               mouse button. Mouse location and axes index values are
               returned

scroll roll - if do_scroll_event is True then these events are returned if
              they occur within an axes. Mouse location and axes index
              values are returned


Dependencies
---------------

matplotlib and wxPython must already be installed.

Example
---------------

An internal example can be seen by calling the module from the command line:

    $ python plot_panel_spectrum.py

Classes
---------------

    PlotPanel() - main call to instantiate a matplotlib canvas in a window
    ZoomSpan() - (internal) sets a mouse button to span zoom in the plot
    ZoomBox() - (internal) sets a mouse button to box zoom in the plot
    CursorSpan() - (internal) sets a mouse button to display a cursor span in the plot
    MiddleEvents() - (internal) toggles events for a middle mouse button if available
    DemoPlotPanel() - full functioning GUI example of the PlotPanel class

History
---------------

This class is an expansion of matplotlib embed in wx example by John Bender
and Edward Abraham, see http://www.scipy.org/Matplotlib_figure_in_a_wx_panel
if it still exists ...

This version allows the user to zoom in on the figure using either a span selector
or a box selector. You can also set a persistent span selector that acts as cursor
references on top of whatever is plotted

ZoomSpan based on matplotlib.widgets.SpanSelector
BoxZoom based on matplotlib.widgets.RectangleSelector
CursorSpan based on matplotlib.widgets.SpanSelector

Author
---------------

Brian J. Soher, Duke University, 20 October, 2010

"""
# Python modules
import math

# 3rd party modules
import matplotlib
import wx

# If we set the backend unconditionally, we sometimes get an undesirable message.
if matplotlib.get_backend() != "WXAgg":
    matplotlib.use('WXAgg')

from matplotlib.transforms import blended_transform_factory
from matplotlib.patches    import Rectangle
from matplotlib.lines      import Line2D




class PlotPanel(wx.Panel):
    """This module contains the PlotPanel class that encapsulates a
    matplotlib canvas within a wxpython widget.

    The PlotPanel has a Figure and a Canvas and 'n' Axes. The user defines
    the  number of axes on Init and this number cannot be changed thereafter.
    However, the user can select which of the n axes are actually displayed
    in the Figure.

    How to interact with the Matplotlib Figure/Axes
    -----------------------------------------------

    This class is a base class for a variety of other more constrained 'plot'
    objects. As such, there is no formal interface for putting data into
    the various Axes, you just keep a reference to the instantiated PlotPanel
    object and access each Axes directly.

    Example (non-functional snippet):

        app   = wx.App( False )
        frame = wx.Frame( None, wx.ID_ANY, 'The PlotPanel', size=(800,800) )
        book  = wx.Notebook(frame, -1, style=wx.BK_BOTTOM)
        panel = wx.Panel(book, -1)
        view = PlotPanel( panel, naxes=2, ... )
        t = np.arange(0.0, 2.0, 0.01)
        sine_data = 1.0*np.sin(2*np.pi*t)
        cos_data = 1.0*np.cos(2*np.pi*t)
        line1 = view.axes[0].plot(t, sine_data, 'b', linewidth=1.0)
        line2 = view.axes[1].plot(t, cos_data, 'g', linewidth=1.0)

    Args:
        parent (wx object): a wx object that the PlotPanel can be inserted into.
            usually a Panel or Sizer.

        naxes (int): default=2, The number of matplotlib axes shown vertically
            in the Figure

        color(tuple, bytes): this is a 3-tuple (r,g,b) or 4-tuple (r,g,b,a) color
            designator that tries to match the matplotlib plot face to the GUI
            color. Default is None, which causes the system color to be retrieved
            and used.

        dpi (int): default=None, Sets matplotlib Figure dpi.

        reversex (bool): default=False, Can be set True to have data plotted
            along reverse direction on X-axis

        zoom (str): default='none', Can be set to 'none', 'span', or 'box'. If
            set to 'span' left button click and drag zooms the figure. A span is
            selected along the x-axis. On release, the axes xlim is adjusted
            accordingly. If zoom mode is 'box', then a zoom box is drawn during
            click and drag and figure is zoomed on both x- and y-axes upon release.
            Clicking in place will cause plot to unzoom to full width.

        reference (bool): default=False, If True, then rignt mouse click and drag
            will draw a span selector in the canvas that persists after release.
            Clicking in place will cause references to reset to non-display.

        middle (bool): default=False, If True, middle mouse events, (press, motion,
            select) will be activated.

        unlink (bool): default=False, If False, Zoom and Reference functionality in
            all plots are syched. Zooming in one will zoom all equally. If True,
            the Zoom and Reference functionality is unique in each plot.

        do_zoom_select_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        do_zoom_motion_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        do_refs_select_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        do_refs_motion_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        do_middle_select_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        do_middle_motion_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        do_middle_press_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        do_scroll_event (bool): default=False, Toggles whether a user defined
            method is called. A default event stub exists on instantiation, but is
            basically a no-op, and should be overloaded by user code.

        uses_collections (bool): default=False. Flag informing if a Collection is
            being used to display data in an Axes. If so, different code is used
            to extract data values during event handling.

        xscale_bump (float): default=0.0, Used to 'nudge' the width of the x-axis
            when resetting zoom. A little bump can be used to ensure that the full
            data plot can be seen left to right.

        yscale_bump (float): default=0.0,  Used to 'nudge' the width of the y-axis
            when resetting zoom. A little bump can be used to ensure that the full
            data plot can be seen top to bottom.

        props_zoom (dict): default=dict(alpha=0.2, facecolor='yellow'), Used to
            set the opacity and color of the zoom span/box as it is being drawn.

        props_cursor (dict): default=dict(alpha=0.2, facecolor='yellow'), Used to
            set the opacity and color of the reference cursor span as it is being
            drawn.

    Returns:
        Nothing, just the object itself

    """

    # Set _EVENT_DEBUG to True to activate printing of messages to stdout during events.
    _EVENT_DEBUG = False

    def __init__(self, parent, naxes=2,
                               color=None, 
                               dpi=None, 
                               reversex=False,
                               zoom='none', 
                               reference=False, 
                               middle=False,
                               unlink=False,
                               do_zoom_select_event=False,
                               do_zoom_motion_event=False,
                               do_refs_select_event=False,
                               do_refs_motion_event=False,
                               do_middle_select_event=False,
                               do_middle_motion_event=False,
                               do_middle_press_event=False,
                               do_scroll_event=False,
                               uses_collections=False,
                               xscale_bump=0.0,
                               yscale_bump=0.0,
                               props_zoom=None,
                               props_cursor=None,
                               **kwargs):

        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.figure import Figure


        # initialize Panel
        if 'id' not in list(kwargs.keys()):
            kwargs['id'] = wx.ID_ANY
        if 'style' not in list(kwargs.keys()):
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )

        self.parent = parent
        self.reversex = reversex
        self.uses_collections = uses_collections        # has to be either Lines OR Collections, not both
        self.unlink = unlink
        self.xscale_bump = xscale_bump
        self.yscale_bump = yscale_bump
        
        # Under GTK track self's size to avoid a continuous flow of size events.
        # For details, see: http://scion.duhs.duke.edu/vespa/project/ticket/28
        self._platform_is_gtk = ("__WXGTK__" in wx.PlatformInfo)
        self._current_size = (-1, -1)

        self.figure = Figure( None, dpi )
        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        #self.figure.set_tight_layout(True)
        
        # Create the N axes, add to figure, also keep a permanent reference to
        # each axes so they can be added or removed from the figure interactively.
        self.axes   = []                    # dynamic list for display
        for i in range(naxes):
            self.axes.append(self.figure.add_subplot(naxes,1,i+1))
        self.naxes = naxes   
        self.all_axes = list(self.axes)     # static list for long term interaction
            
        self.zoom = []
        self.refs = []
        self.middle = []
        
        self.set_color( color )
        self._set_size()
        self._resizeflag = False

        self.Bind(wx.EVT_IDLE, self._on_idle)
        self.Bind(wx.EVT_SIZE, self._on_size)

        # ensure that default properties exist
        if not props_zoom:
            props_zoom = dict(alpha=0.2, facecolor='yellow')

        if not props_cursor:
            props_cursor = dict(alpha=0.2, facecolor='purple')

        #----------------------------------------------------------------------
        # enable Zoom, Reference, Middle and Scroll functionality as required

        if zoom == 'span':
        
            if not unlink:
                self.zoom = ZoomSpan( self, self.all_axes,
                                      useblit=True,
                                      do_zoom_select_event=do_zoom_select_event,
                                      do_zoom_motion_event=do_zoom_motion_event,
                                      rectprops=props_zoom)
            else:
                for axes in self.axes:
                    self.zoom.append( ZoomSpan( self, [axes],
                                          useblit=True,
                                          do_zoom_select_event=do_zoom_select_event,
                                          do_zoom_motion_event=do_zoom_motion_event,
                                          rectprops=props_zoom))
        if zoom == 'box':
            
            if not unlink:
                self.zoom = ZoomBox(  self, self.axes,
                                      drawtype='box',
                                      useblit=True,
                                      button=1,
                                      do_zoom_select_event=do_zoom_select_event,
                                      do_zoom_motion_event=do_zoom_motion_event,
                                      spancoords='data',
                                      rectprops=props_zoom)
            else:
                for axes in self.axes:
                    self.zoom.append(ZoomBox(  self, [axes],
                                          drawtype='box',
                                          useblit=True,
                                          button=1,
                                          do_zoom_select_event=do_zoom_select_event,
                                          do_zoom_motion_event=do_zoom_motion_event,
                                          spancoords='data',
                                          rectprops=props_zoom))
        if reference:
            
            if not unlink:
                self.refs = CursorSpan(self, self.axes,
                                      useblit=True,
                                      do_refs_select_event=do_refs_select_event,
                                      do_refs_motion_event=do_refs_motion_event,
                                      rectprops=props_cursor)
            else:
                for axes in self.axes:
                    self.refs.append(CursorSpan(self, [axes],
                                          useblit=True,
                                          do_refs_select_event=do_refs_select_event,
                                          do_refs_motion_event=do_refs_motion_event,
                                          rectprops=props_cursor))
        if middle:

            if not unlink:
                self.middle = MiddleEvents(self, self.axes,
                                      do_middle_select_event=do_middle_select_event,
                                      do_middle_motion_event=do_middle_motion_event,
                                      do_middle_press_event=do_middle_press_event)
            else:
                for axes in self.axes:
                    self.middle.append( MiddleEvents(self, [axes],
                                          do_middle_select_event=do_middle_select_event,
                                          do_middle_motion_event=do_middle_motion_event))

        self.do_motion_event = True  
        self.motion_id = self.canvas.mpl_connect('motion_notify_event', self._on_move)
        
        self.do_scroll_event = do_scroll_event
        if self.do_scroll_event:
            self.scroll_id = self.canvas.mpl_connect('scroll_event', self._on_scroll)



    #=======================================================
    #
    #           Internal Helper Functions  
    #
    #=======================================================

    def _dprint(self, a_string):
        if self._EVENT_DEBUG:
            print(a_string)

    def _on_size( self, event ):
        if self._platform_is_gtk:
            # This is a workaround for linux bug where 1 pixel sized change creates loop
            current_x, current_y = self._current_size
            new_x, new_y = tuple(event.GetSize())
            if (abs(current_x - new_x) > 1) or (abs(current_y - new_y) > 1):
                self._resizeflag = True
            else:
                event.Skip(False)       # Size change <= 1 pixel, ignore it.
        else:
            self._resizeflag = True

    def _on_idle( self, evt ):
        if self._resizeflag:
            self._resizeflag = False
            self._set_size()

    def _set_size( self ):
        pixels = tuple( self.parent.GetClientSize() )
        self.SetSize( pixels )
        self.canvas.SetSize( pixels )
        self.figure.set_size_inches( float( pixels[0] )/self.figure.get_dpi(),
                                     float( pixels[1] )/self.figure.get_dpi() )
        self._current_size = pixels

    def _on_move(self, event):
        """
        Internal method organizes data sent to an external user defined event
        handler for motion events. We gather values from either Collection or
        Line plots, determine which axis we are in, then call on_motion()

        """
        if event.inaxes == None or not self.do_motion_event: return
        
        bounds = event.inaxes.dataLim.bounds
        value = self.get_values(event)
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i

        self.on_motion(event.xdata, event.ydata, value, bounds, iaxis)
        
    def _on_scroll(self, event):
        """
        Internal method organizes data sent to an external user defined event
        handler for scroll events. We determine which axis we are in, then
        call on_scroll()

        """
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i
                
        self.on_scroll(event.button, event.step, iaxis)        


    def get_values(self, event):
        """
        Included in plot_panel object so users can overwrite it if necessary.
        Default functionality: Determine which axes the mouse is in, return a
        list of data values at the x location of the cursor.
        
        This is complicated by the fact that some plot_panels put their plots
        into the lines attribute, while others use the collections attribute.
        This is solved by designating this distinction using a flag at the 
        main level. Caveat, plot_panels can have only one or the other types 
        of plots, line or collection based.
        
        """
        value = []
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds

        if not all([math.isfinite(x) for x in [x0, y0, x1, y1]]):    # typical when plot is empty
            return [0.0,]
    
        # sporadic event calls yield 'event.inaxes == None' so we test for that here
        if event.inaxes:
            if self.uses_collections:
                
                # Expects LineCollection - Figures out the offset of the cursor and
                # returns the path value for the x-position of that line plot.
                # Note. Not guaranteed to work for any other types of collections
                if event.inaxes.collections:
                    npts = len(event.inaxes.collections[0]._paths[0].vertices[:,1])
                    indx = int(round((npts-1) * (event.xdata-x0)/x1))
                    if self.reversex:   indx = npts - indx - 1
                    if indx > (npts-1): indx = npts - 1
                    if indx < 0:        indx = 0
                    for i,path in enumerate(event.inaxes.collections[0]._paths):
                        offset = event.inaxes.collections[0]._uniform_offsets[i,1]
                        dat = path.vertices[:,1]
                        if indx < len(dat):
                            value.append(dat[indx]-offset)
            else:
                if event.inaxes.lines:
                    data = event.inaxes.lines[0].get_ydata()
                    if len(data.shape)==0:
                        value = [data]
                    else:
                        npts = len(event.inaxes.lines[0].get_ydata())
                        indx = int(round((npts-1) * (event.xdata-x0)/x1))
                        if self.reversex:   indx = npts - indx - 1
                        if indx > (npts-1): indx = npts - 1
                        if indx < 0:        indx = 0
                        for line in event.inaxes.lines:
                            dat = line.get_ydata()
                            if indx < len(dat):
                                value.append(dat[indx])
        if value == []: 
            value = [0.0,]
            
        return value


    #=======================================================
    #
    #           User Accessible Functions  
    #
    #=======================================================

    def set_color( self, rgbtuple=None ):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()
        clr = [c/255. for c in rgbtuple]
        self.figure.set_facecolor( clr )
        self.figure.set_edgecolor( clr )
        self.canvas.SetBackgroundColour( wx.Colour( *rgbtuple ) )

    def refresh_cursors(self):
        """ redraws the reference cursor span on user request """
        if self.refs == None: return
        if self.refs.rect == []: return
        for axes, rect in zip(self.axes, self.refs.rect):
            axes.add_patch(rect)
        self.canvas.draw()

    def change_naxes(self, n):
        """
        Allow user to determine interactively which of the N axes are included
        in the figure. NB. This method irrevocably removes axes from the figure.
        User supplies only the number of axes to include and the first 1:n
        axes in the long term storage list are retaine in the figure. This method
        also updates the axes lists in any zoom, refs or middle objects.
        
        """
        if n > self.naxes or n < 0 or n == len(self.figure.axes):
            return
        
        self.axes = self.all_axes[0:n]  

        figure_axes = list(self.figure.axes)
        for axes in figure_axes:
            self.figure.delaxes(axes)
            
        for axes in self.axes:
            ax = self.figure.add_axes(axes)

        if not self.unlink:
            if self.zoom:
                self.zoom.axes = self.axes
                 
            if self.refs:
                self.refs.axes = self.axes
     
            if self.middle:
                self.middle.axes = self.axes

        self.canvas.draw()

    def display_naxes(self, flags):
        """
        Specifiy which of the N axes in the figure are drawn. 'Flags' is a
        boolean list the same length as the self.all_axes list. Axes that
        correspond to flags set to True are drawn. NB. No axes are destroyed
        by this method. This method also updates the axes lists in any zoom,
        refs or middle objects.
        
        """
        ncurrent = len(self.all_axes)
        nflags = len(flags)
        if nflags != ncurrent: return

        faxes = list(self.figure.axes)
        for axes in faxes:
            self.figure.delaxes(axes)
            
        for i, axes in enumerate(self.axes):
            if flags[i] != False:
                ax = self.figure.add_axes(axes)

        if not self.unlink:
            if self.zoom:
                self.zoom.axes = self.axes
                 
            if self.refs:
                self.refs.axes = self.axes
     
            if self.middle:
                self.middle.axes = self.axes

        self.canvas.draw()

    def new_axes(self, axes):
        if isinstance(axes, list):
            self.axes = axes
        elif isinstance(axes, matplotlib.axes.Axes):
            self.axes = [axes]
        else:
            return

        if self.zoom is not None:
            self.zoom.new_axes(self.axes)
        if self.reference is not None:
            self.refs.new_axes(self.axes)
            
        if self.canvas is not self.axes[0].figure.canvas:
            self.canvas.mpl_disconnect(self.motion_id)
            self.canvas = self.axes[0].figure.canvas
            self.motion_id = self.canvas.mpl_connect('motion_notify_event', self._on_move)       
        if self.figure is not self.axes[0].figure:
            self.figure = self.axes[0].figure

        
    #=======================================================
    #
    #           Default Event Handlers  
    #
    #=======================================================
        
    def on_motion(self, xdata, ydata, value, bounds, iaxis):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_move, xdata='+str(xdata)+'  ydata='+str(ydata)+'  val='+str(value)+'  bounds = '+str(bounds)+'  iaxis='+str(iaxis))
        
    def on_scroll(self, button, step, iaxis):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_move, button='+str(button)+'  step='+str(step)+'  iaxis='+str(iaxis))
        
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_zoom_select, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val)+'  ymin='+str(ymin)+'  ymax='+str(ymax))
        
    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_zoom_move, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val)+'  ymin='+str(ymin)+'  ymax='+str(ymax))

    def on_refs_select(self, xmin, xmax, val, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_refs_select, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val))
        
    def on_refs_motion(self, xmin, xmax, val, iplot=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_refs_move, xmin='+str(xmin)+'  xmax='+str(xmax)+'  val='+str(val))

    def on_middle_select(self, xstr, ystr, xend, yend, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('ext on_middle_select, X(str,end)='+str(xstr)+','+str(xend)+'  Y(str,end)='+str(ystr)+','+str(yend)+'  Index = '+str(indx))
        
    def on_middle_motion(self, xcur, ycur, xprev, yprev, indx):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_middle_move, X(cur,prev)='+str(xstr)+','+str(xprev)+'  Y(cur,prev)='+str(ystr)+','+str(yprev)+'  Index = '+str(indx))

    def on_middle_press(self, xloc, yloc, indx, bounds=None, xdata=None, ydata=None):
        """ placeholder, overload for user defined event handling """
        self._dprint('on_middle_press, Xloc='+str(xloc)+'  Yloc='+str(yloc)+'  Index = '+str(indx))
        

class ZoomSpan:
    """ Select a min/max range of the x axis for a matplotlib Axes """

    def __init__(self, parent, axes, 
                               minspan=None, 
                               useblit=False, 
                               rectprops=None, 
                               do_zoom_select_event=False, 
                               do_zoom_motion_event=False):
        """
        Create a span selector in axes.  When a selection is made, clear
        the span and call onselect with

          onselect(vmin, vmax)

        If minspan is not None, ignore events smaller than minspan

        The span rect is drawn with rectprops; default
          rectprops = dict(facecolor='yellow', alpha=0.2)

        set the visible attribute to False if you want to turn off
        the functionality of the span selector

        """
        if rectprops is None:
            rectprops = dict(facecolor='yellow', alpha=0.2)

        self.parent = parent
        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []

        self.rect = []
        self.background = None
        self.pressv = None

        self.rectprops = rectprops
        self.do_zoom_select_event = do_zoom_select_event
        self.do_zoom_motion_event = do_zoom_motion_event
        self.useblit = useblit
        self.minspan = minspan

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prev = (0, 0)

        self.new_axes(axes)


    def new_axes(self,axes):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas

            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))
        
        for axes in self.axes:
            trans = blended_transform_factory(axes.transData, axes.transAxes)
            self.rect.append(Rectangle( (0,0), 0, 1,
                                   transform=trans,
                                   visible=False,
                                   **self.rectprops ))
        if not self.useblit:
            for axes, rect in zip(self.axes, self.rect):
                axes.add_patch(rect)

    def update_background(self, event):
        '''force an update of the background'''
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        '''return True if event should be ignored'''
        return  event.inaxes not in self.axes or not self.visible or event.button !=1

    def press(self, event):
        '''on button press event'''
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = False

        # remove the dynamic artist(s) from background bbox(s)
        for axes, rect in zip(self.axes, self.rect):
            if rect in axes.patches:
                axes.patches.remove(rect)
                self.canvas.draw()

        for rect in self.rect:
            rect.set_visible(self.visible)
            
        self.pressv = event.xdata
        return False

    def release(self, event):
        '''on button release event'''
        if self.pressv is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = True

        for rect in self.rect:
            rect.set_visible(False)

        # left-click in place resets the x-axis 
        if event.xdata == self.pressv:
            x0, y0, x1, y1 = event.inaxes.dataLim.bounds  
            if all([math.isfinite(x) for x in [x0, y0, x1, y1]]):  # typical when plot is empty
                xdel = self.parent.xscale_bump * (x1 - x0)
                ydel = self.parent.yscale_bump * (y1 - y0)
                for axes in self.axes:
                    if self.parent.reversex:
                        axes.set_xlim(x0+x1+xdel,x0-xdel)
                    else:
                        axes.set_xlim(x0-xdel,x0+x1+xdel)
                    axes.set_ylim(y0-ydel,y0+y1+ydel)
                self.canvas.draw()

                if self.do_zoom_select_event:
                    self.parent.on_zoom_select(x0-xdel, x0+x1+xdel, [0.0], y0-ydel, y0+y1+ydel, reset=True)
            
            return

        vmin = self.pressv
        vmax = event.xdata or self.prev[0]

        if vmin>vmax: vmin, vmax = vmax, vmin
        span = vmax - vmin
        if self.minspan is not None and span<self.minspan: return

        for axes in self.axes:
            if self.parent.reversex:
                axes.set_xlim((vmax, vmin))
            else:
                axes.set_xlim((vmin, vmax))

        self.canvas.draw()

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes is not None:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_zoom_select_event and data_test:
            # gather the values to report in a selection event
            value = self.parent.get_values(event)
            self.parent.on_zoom_select(vmin, vmax, value, None, None) 
            
        self.pressv = None
        
        return False

    def update(self):
        '''draw using newfangled blit or oldfangled draw depending on useblit'''
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, rect in zip(self.axes, self.rect):
                axes.draw_artist(rect)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()

        return False

    def onmove(self, event):
        '''on motion notify event'''
        if self.pressv is None or self.ignore(event): return
        x, y = event.xdata, event.ydata
        self.prev = x, y

        minv, maxv = x, self.pressv
        if minv>maxv: minv, maxv = maxv, minv
        for rect in self.rect:
            rect.set_x(minv)
            rect.set_width(maxv-minv)

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_zoom_motion_event and data_test: 
            vmin = self.pressv
            vmax = event.xdata or self.prev[0]
            if vmin>vmax: vmin, vmax = vmax, vmin
            value = self.parent.get_values(event)
            self.parent.on_zoom_motion(vmin, vmax, value, None, None) 

        self.update()
        return False




class CursorSpan:
    """ Indicate two vertical reference lines along a matplotlib Axes """

    def __init__(self, parent, axes,
                               minspan=None, 
                               useblit=False, 
                               rectprops=None, 
                               do_refs_select_event=False, 
                               do_refs_motion_event=False):
        """
        Create a span selector in axes. When a selection is made, clear the
        span and call onselect with

          onselect(vmin, vmax)

        and clear the span.

        If minspan is not None, ignore events smaller than minspan

        The span rect is drawn with rectprops; default
          rectprops = dict(facecolor='none', alpha=1.0)

        set the visible attribute to False if you want to turn off
        the functionality of the span selector

        """
        if rectprops is None:
            rectprops = dict(facecolor='none')

        self.parent = parent
        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []

        self.rect = []
        self.background = None
        self.pressv = None

        self.rectprops = rectprops
        self.do_refs_select_event = do_refs_select_event
        self.do_refs_motion_event = do_refs_motion_event
        self.useblit = useblit
        self.minspan = minspan

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prev  = (0,0)

        self.new_axes(axes)


    def set_span(self, xmin, xmax):
        
        x0, y0, x1, y1 = self.axes[0].dataLim.bounds
        if all([math.isfinite(x) for x in [x0, y0, x1, y1]]):  # typical when plot is empty
            self.visible = True
            if xmin < x0: xmin = x0
            if xmax > (x0+x1): xmax = x0+x1

            for rect in self.rect:
                rect.set_x(xmin)
                rect.set_width(xmax-xmin)

            self.canvas.draw()
        

    def new_axes(self,axes):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas

            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))
        
        for axes in self.axes:
            trans = blended_transform_factory(axes.transData, axes.transAxes)
            self.rect.append(Rectangle( (0,0), 0, 1,
                                   transform=trans,
                                   visible=False,
                                   **self.rectprops ))

        if not self.useblit: 
            for axes, rect in zip(self.axes, self.rect):
                axes.add_patch(rect)


    def update_background(self, event):
        '''force an update of the background'''
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        '''return True if event should be ignored'''
        return  event.inaxes not in self.axes or not self.visible or event.button !=3

    def press(self, event):
        '''on button press event'''
        self.visible = True
        if self.ignore(event): return
        self.buttonDown = True
        
        # only send one motion event while selecting
        if self.do_refs_motion_event:
            self.parent.do_motion_event = False

        # remove the dynamic artist(s) from background bbox(s)
        for axes, rect in zip(self.axes, self.rect):
            if rect in axes.patches:
                axes.patches.remove(rect)
                self.canvas.draw()
        
        for rect in self.rect:
            rect.set_visible(self.visible)
        self.pressv = event.xdata
        return False

    def release(self, event):
        '''on button release event'''
        if self.pressv is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_refs_motion_event:
            self.parent.do_motion_event = True

        # left-click in place resets the x-axis 
        if event.xdata == self.pressv:
            self.visible = not self.visible
            for axes, rect in zip(self.axes, self.rect):
                rect.set_visible(self.visible)
                axes.add_patch(rect)
            self.canvas.draw()  
            self.pressv = None          
            return

        vmin = self.pressv
        vmax = event.xdata or self.prev[0]

        if vmin>vmax: vmin, vmax = vmax, vmin
        span = vmax - vmin
        # don't add reference span, if min span not achieved
        if self.minspan is not None and span<self.minspan: return

        for axes, rect in zip(self.axes, self.rect):
            rect.set_visible(True)
            axes.add_patch(rect)
            self.canvas.draw()

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_refs_select_event and data_test:
            # don't gather values if no onselect event
            value = self.parent.get_values(event)
            self.parent.on_refs_select(vmin, vmax, value)
        
        self.pressv = None
        
        return False

    def update(self):
        '''draw using newfangled blit or oldfangled draw depending on useblit'''
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, rect in zip(self.axes, self.rect):
                axes.draw_artist(rect)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()

        return False

    def onmove(self, event):
        '''on motion notify event'''
        if self.pressv is None or self.ignore(event): return
        x, y = event.xdata, event.ydata
        self.prev = x, y

        minv, maxv = x, self.pressv
        if minv>maxv: minv, maxv = maxv, minv
        for rect in self.rect:
            rect.set_x(minv)
            rect.set_width(maxv-minv)

        # sporadic event calls yield event.inaxes == None so we test here
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_refs_motion_event and data_test:
            vmin = self.pressv
            vmax = event.xdata or self.prev[0]
            if vmin>vmax: vmin, vmax = vmax, vmin
            value = self.parent.get_values(event)            
            self.parent.on_refs_motion(vmin, vmax, value) 

        self.update()
        return False




class ZoomBox:
    """ Select a min/max range of the x and y axes for a matplotlib Axes """

    def __init__(self, parent, axes,
                             drawtype='box',
                             minspanx=None, 
                             minspany=None, 
                             useblit=False,
                             lineprops=None, 
                             rectprops=None,
                             do_zoom_select_event=False, 
                             do_zoom_motion_event=False,
                             spancoords='data',
                             button=None):

        """
        Create a selector in axes.  When a selection is made, clear
        the span and call onselect with

          onselect(pos_1, pos_2)

        and clear the drawn box/line. There pos_i are arrays of length 2
        containing the x- and y-coordinate.

        If minspanx is not None then events smaller than minspanx
        in x direction are ignored(it's the same for y).

        The rect is drawn with rectprops; default
          rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.5, fill=False)

        The line is drawn with lineprops; default
          lineprops = dict(color='black', linestyle='-',
                           linewidth = 2, alpha=0.5)

        Use type if you want the mouse to draw a line, a box or nothing
        between click and actual position ny setting

        drawtype = 'line', drawtype='box' or drawtype = 'none'.

        spancoords is one of 'data' or 'pixels'.  If 'data', minspanx
        and minspanx will be interpreted in the same coordinates as
        the x and y axis, if 'pixels', they are in pixels

        button is a list of integers indicating which mouse buttons should
        be used for rectangle selection.  You can also specify a single
        integer if only a single button is desired.  Default is None, which
        does not limit which button can be used.
        Note, typically:
         1 = left mouse button
         2 = center mouse button (scroll wheel)
         3 = right mouse button

        """
        self.parent = parent
        self.axes = None
        self.canvas = None
        self.visible = True
        self.cids = []
        
        self.active = True                    # for activation / deactivation
        self.to_draw = []
        self.background = None

        self.do_zoom_select_event = do_zoom_select_event
        self.do_zoom_motion_event = do_zoom_motion_event
        
        self.useblit = useblit
        self.minspanx = minspanx
        self.minspany = minspany

        if button is None or isinstance(button, list):
            self.validButtons = button
        elif isinstance(button, int):
            self.validButtons = [button]

        assert(spancoords in ('data', 'pixels'))

        self.spancoords = spancoords
        self.eventpress = None          # will save the data (position at mouseclick)
        self.eventrelease = None        # will save the data (pos. at mouserelease)

        self.new_axes(axes, rectprops)


    def new_axes(self,axes, rectprops=None):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas
            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
            self.cids.append(self.canvas.mpl_connect('draw_event', self.update_background))
            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            
        if rectprops is None:
            rectprops = dict(facecolor='white', 
                             edgecolor= 'black',
                             alpha=0.5, 
                             fill=False)
        self.rectprops = rectprops

        for axes in self.axes:
            self.to_draw.append(Rectangle((0,0), 0, 1,visible=False,**self.rectprops))

        for axes,to_draw in zip(self.axes, self.to_draw):
            axes.add_patch(to_draw)
    

    def update_background(self, event):
        '''force an update of the background'''
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def ignore(self, event):
        '''return True if event should be ignored'''
        # If ZoomBox is not active :
        if not self.active:
            return True

        # If canvas was locked
        if not self.canvas.widgetlock.available(self):
            return True

        # Only do selection if event was triggered with a desired button
        if self.validButtons is not None:
            if not event.button in self.validButtons:
                return True

        # If no button pressed yet or if it was out of the axes, ignore
        if self.eventpress == None:
            return event.inaxes not in self.axes

        # If a button pressed, check if the release-button is the same
        return  (event.inaxes not in self.axes or
                 event.button != self.eventpress.button)

    def press(self, event):
        '''on button press event'''
        # Is the correct button pressed within the correct axes?
        if self.ignore(event): return
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = False


        # make the drawn box/line visible get the click-coordinates, button, ...
        for to_draw in self.to_draw:
            to_draw.set_visible(self.visible)
        self.eventpress = event
        return False


    def release(self, event):
        '''on button release event'''
        if self.eventpress is None or self.ignore(event): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas
        
        # only send one motion event while selecting
        if self.do_zoom_motion_event:
            self.parent.do_motion_event = True
        
        # make the box/line invisible again
        for to_draw in self.to_draw:
            to_draw.set_visible(False)
        
        # left-click in place resets the x-axis or y-axis
        if self.eventpress.xdata == event.xdata and self.eventpress.ydata == event.ydata:
            x0, y0, x1, y1 = event.inaxes.dataLim.bounds
            if all([math.isfinite(x) for x in [x0,y0,x1,y1]]):  # typical when plot is empty
                xdel = self.parent.xscale_bump * (x1 - x0)
                ydel = self.parent.yscale_bump * (y1 - y0)
                for axes in self.axes:
                    if self.parent.reversex:
                        axes.set_xlim(x0+x1+xdel,x0-xdel)
                    else:
                        axes.set_xlim(x0-xdel,x0+x1+xdel)
                    axes.set_ylim(y0-ydel,y0+y1+ydel)
                self.canvas.draw()

                if self.do_zoom_select_event:
                    self.parent.on_zoom_select(x0-xdel, x0+x1+xdel, [0.0], y0-ydel, y0+y1+ydel, reset=True)
            
            return
        
        self.canvas.draw()
        self.eventrelease = event   # release coordinates, button, ...

        if self.spancoords=='data':
            xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
            xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata
            # calculate dimensions of box or line get values in the right order
        elif self.spancoords=='pixels':
            xmin, ymin = self.eventpress.x, self.eventpress.y
            xmax, ymax = self.eventrelease.x, self.eventrelease.y
        else:
            raise ValueError('spancoords must be "data" or "pixels"')

        # assure that min<max values and x and y values are not equal
        if xmin>xmax: xmin, xmax = xmax, xmin
        if ymin>ymax: ymin, ymax = ymax, ymin
        if xmin == xmax: xmax = xmin*1.0001
        if ymin == ymax: ymax = ymin*1.0001

        spanx = xmax - xmin
        spany = ymax - ymin
        xproblems = self.minspanx is not None and spanx<self.minspanx
        yproblems = self.minspany is not None and spany<self.minspany
        if (xproblems or  yproblems):
            """Box too small"""    # check if drawed distance (if it exists) is
            return                 # not to small in neither x nor y-direction
        
        for axes in self.axes:
            if self.parent.reversex:
                axes.set_xlim((xmax,xmin))
            else:
                axes.set_xlim((xmin,xmax))
            axes.set_ylim((ymin,ymax))
        self.canvas.draw()

        # sporadic event calls yield event.inaxes == None so we test here
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]

        if self.do_zoom_select_event and data_test:
            value = self.parent.get_values(event)
            self.parent.on_zoom_select(xmin, xmax, value, ymin, ymax) # zeros are for consistency with box zoom

        self.eventpress   = None              # reset variables to inital values
        self.eventrelease = None
        
        return False


    def update(self):
        '''draw using newfangled blit or oldfangled draw depending on useblit'''
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for axes, to_draw in zip(self.axes, self.to_draw):
                axes.draw_artist(to_draw)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()
        return False


    def onmove(self, event):
        '''on motion notify event if box/line is wanted'''
        if self.eventpress is None or self.ignore(event): return
        x,y = event.xdata, event.ydata              # actual position (with
                                                    #   (button still pressed)
        minx, maxx = self.eventpress.xdata, x       # click-x and actual mouse-x
        miny, maxy = self.eventpress.ydata, y       # click-y and actual mouse-y
        if minx>maxx: minx, maxx = maxx, minx       # get them in the right order
        if miny>maxy: miny, maxy = maxy, miny       
        for to_draw in self.to_draw:
            to_draw.set_x(minx)                    # set lower left of box
            to_draw.set_y(miny)
            to_draw.set_width(maxx-minx)           # set width and height of box
            to_draw.set_height(maxy-miny)

        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        data_test = False
        if event.inaxes:
            if self.parent.uses_collections:
                data_test = event.inaxes.collections!=[]
            else:
                data_test = event.inaxes.lines!=[]
        
        if self.do_zoom_motion_event and data_test:
            # gather the values to report in a selection event
            value = self.parent.get_values(event)
            self.parent.on_zoom_motion(minx, maxx, value, miny, maxy) # zeros are for consistency with box zoom
        
        self.update()
        return False

    def set_active(self, active):
        ''' Use to de/activate RectangleSelector with a boolean variable 'active' '''
        self.active = active

    def get_active(self):
        ''' to get status of active mode (boolean variable)'''
        return self.active



class MiddleEvents:
    """ Report events having to do with the middle button to a user method """

    def __init__(self, parent, axes, 
                               do_middle_select_event=False, 
                               do_middle_motion_event=False,
                               do_middle_press_event=False):
        """
        Events can be set up granularly to call some or all of 'select',
        'motion', and/or 'press' events.

        If you set these events to be On, then you are responsible for
        overloading the appropriate method to deal with them.

        """
        self.parent = parent
        self.axes = None
        self.canvas = None
        self.cids = []

        self.background = None
        self.pressxy = None
        self.axes_index = None

        self.do_middle_select_event = do_middle_select_event
        self.do_middle_motion_event = do_middle_motion_event
        self.do_middle_press_event  = do_middle_press_event

        # Needed when dragging out of axes
        self.buttonDown = False
        self.prevxy = (0,0)

        self.new_axes(axes)


    def new_axes(self,axes):
        self.axes = axes
        if self.canvas is not axes[0].figure.canvas:
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)

            self.canvas = axes[0].figure.canvas

            self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))
            self.cids.append(self.canvas.mpl_connect('button_press_event', self.press))
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.release))
    
    def ignore(self, event):
        '''return True if event should be ignored'''
        return  event.inaxes not in self.axes or event.button !=2

    def press(self, event):
        '''on button press event'''
        if self.ignore(event): return
        self.buttonDown = True
        
        for i in range(len(self.axes)):
            if event.inaxes == self.axes[i]:
                self.axes_index = i
        
        # only send one motion event while selecting
        if self.do_middle_motion_event:
            self.parent.do_motion_event = False

        bounds = event.inaxes.dataLim.bounds

        self.pressxy = event.x, event.y
        self.prevxy  = event.x, event.y

        if self.do_middle_press_event:
            self.parent.on_middle_press(event.x, event.y, self.axes_index, bounds=bounds, xdata=event.xdata, ydata=event.ydata)
        
        return False

    def release(self, event):
        '''on button release event'''
        if self.pressxy is None or (self.ignore(event) and not self.buttonDown): return

        self.parent.SetFocus()  # sets focus into Plot_Panel widget canvas

        self.buttonDown = False
        
        # only send one motion event while selecting
        if self.do_middle_motion_event:
            self.parent.do_motion_event = True

        xstr, ystr = self.pressxy
        xend = event.x
        yend = event.y

        if self.do_middle_select_event:
            self.parent.on_middle_select(xstr, ystr, xend, yend, self.axes_index)

        self.axes_index = None
        self.pressxy = None
        return False

    def onmove(self, event):
        '''on motion notify event'''
        if self.pressxy is None or self.ignore(event): return
        xcurrent, ycurrent = event.x, event.y
        xprevious, yprevious = self.prevxy
        
        self.prevxy = event.x, event.y

        if self.do_middle_motion_event:
            self.parent.on_middle_motion(xcurrent, ycurrent, xprevious, yprevious, self.axes_index) 

        return False
    


#------------------------------------------------
# Test Code
#------------------------------------------------

if __name__ == '__main__':

    import numpy as np

    class DemoPlotPanel(PlotPanel):
        """Plots several lines in distinct colors."""

        # Activate event messages
        _EVENT_DEBUG = True

        def __init__( self, parent, **kwargs ):
            # initiate plotter
            PlotPanel.__init__( self, parent, **kwargs )  
            self.parent = parent


    app   = wx.App( False )
    frame = wx.Frame( None, wx.ID_ANY, 'WxPython and Matplotlib - PlotPanel', size=(800,800) )
    
    nb = wx.Notebook(frame, -1, style=wx.BK_BOTTOM)
    
    panel1 = wx.Panel(nb, -1)
    view = DemoPlotPanel( panel1, naxes=2,
                                  color = (122,122,122,255),
                                  zoom='span', 
                                  reference=True,
                                  middle=True,
                                  unlink=False,
                                  do_zoom_select_event=True,
                                  do_zoom_motion_event=True,
                                  do_refs_select_event=True,
                                  do_refs_motion_event=True,
                                  do_middle_select_event=True,
                                  do_middle_motion_event=True,
                                  do_scroll_event=True )

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
    panel1.SetSizer(sizer)
    view.Fit()    

    nb.AddPage(panel1, "One")
    
    t = np.arange(0.0, 2.0, 0.01)
    s = 1.0*np.sin(2*np.pi*t)
    c = 1.0*np.cos(2*np.pi*t)

    view.set_color( (255,255,255) )

    line1 = view.axes[0].plot(t, s, 'b', linewidth=1.0)
    line2 = view.axes[1].plot(t, c, 'g', linewidth=3.0)
    view.axes[0].set_xlabel('time (s)')
    for axes in view.axes:
        axes.set_ylabel('voltage (mV)')
        axes.grid(True)
    view.axes[0].set_title('About as simple as it gets, folks')
    
    view.canvas.figure.subplots_adjust(left=0.20,
                                         right=0.95,
                                         bottom=0.15,
                                         top=0.95,
                                         wspace=0.0,
                                         hspace=0.01)
    view.canvas.draw()
    frame.Show()
    app.MainLoop()


