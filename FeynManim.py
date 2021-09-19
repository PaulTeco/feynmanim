from manim import *


class FeynManim:

    """
    A class to construct parts (curves) of Feynman Diagrams with Manim


    Main Attributes
    ----------
    geometry : lambda function
        parametrizes the geometry of the curve (default is a straight line)
    line_style : lambda function
        parametrizes the style of the curve (default is plain)
    flow_tube : Manim StreamLines Object
        can be used to display flow
    mobject : Manim ParametricCurve or VGroup Object
        can be used to display the curve with Manim

    Geometric Attributes
    ----------
    vertex1 : numpy.array([float, float, float])
        coordinates of the start of the curve
    vertex2 : numpy.array([float, float, float])
        coordinates of the start of the curve
    geolength : float
        length of the geometry of the curve
    center : numpy.array([float, float, float])
        coordinates of the center of the curve if it is an arc
    radius : float
        radius of the curve if it is an arc

    Additional Attributes
    ----------
    styled_geometry : lambda function
        a combination of the geometry and the style function
    momentum_flow_field : lambda function
        a realization of the momentum flow as a vector field

    Methods
    -------
    set_geometry(geotype,invert=False,arclen=0.5)
        sets the geometry of the curve
    set_style(style,wiggles=False,amp=1)
        sets the style of the curve
    add_arrow(flip=False,pos=0.5,color=WHITE)
        adds a triangle on a curve to represent e.g. Fermion line arrows
    build_flow_tube(color=WHITE,delta_x=0.3,delta_y=0.3,stroke_width=2,max_anchors_per_line=30,cut_start_bleed=0,cut_end_bleed=0,invert=False)
        uses manim's StreamLines method to create a flowing animation on a curve
    build_mobject(color=WHITE)
        builds the object as an instance of manim's ParametricFunction or VGroup
    """

    def __init__(self,p1,p2):

        """Initiates the Object as a plain line

        Parameters
        ----------
        p1 : numpy.array([float, float, float])
            coordinates of the start of the curve
        p2 : numpy.array([float, float, float])
            coordinates of the end of the curve
        """

        self.vertex1 = np.array(p1)
        self.vertex2 = np.array(p2)
        self.arrows = []

        # initiate a simple line without decorators
        self.set_geometry("line")
        self.set_style("plain")


    # # Helping Functions
    def _derivative(self,func):
        """returns a derivative of a lambda function as another lambda function"""

        delta_t = 0.0000000001
        return lambda t: (func(t+delta_t)-func(t-delta_t))/(2*delta_t)

    def _orthonormal(self,func):
        """for a R^1 -> R^3 lambda function, returns the orthonormal projection in the xy-plane"""

        return lambda t: np.array([func(t)[1],-func(t)[0],0])/np.linalg.norm(func(t))


    def _orthonormal_vec(self,v):
        """for a vector of R^3, returns the orthonormal projection in the xy-plane"""
        return np.array([-v[1],v[0],0])

    def _masking_function(self,pos,geo,thickness=0.5,curve_pts=50,cut_start_bleed=0,cut_end_bleed=0):
        """
        takes a geometry and position and evaluates as 1 close (thickness) to the geometry and 0 elsewhere
        this can be used to modify a vector field so it is 0 everywhere not close enough to the function
        cut_start_bleed and cut_end_bleed take percentages from the start/end of the flow tube (0-1)
        """
        
        curve_dots = [geo( (t*(1-cut_start_bleed-cut_end_bleed)+cut_start_bleed) / (curve_pts-1) ) for t in range(curve_pts)]

        if(np.min([np.linalg.norm(v-pos) for v in curve_dots]) <= thickness):
            return 1
        else: 
            return 0


    # find approximate xmin xmax ymin ymax of the curve. this considerably reduces computing time for the flow field
    def _ext_of_curve(self,curve,param_range=[0,1],pts=100):
        """for a given curve and parameter range, approximates xmin xmax ymin ymax of the curve by evaluating pts evenly spaced points"""

        t_pts = np.linspace(*param_range,pts)
        sample = np.array([curve(t) for t in t_pts])
        xmin, xmax = [np.min(sample[:,0]), np.max(sample[:,0])]
        ymin, ymax = [np.min(sample[:,1]), np.max(sample[:,1])]
        return [xmin,xmax,ymin,ymax]


    # # # set geometry of the curve ------------------------------------------
    def set_geometry(self,geotype,invert=False,arclen=0.5): 

        """Set the geometry of the curve

        Parameters
        ----------
        geotype : str
            either "arc" or "line"
        invert : bool
            A flag used to switch start and end; mirror arcs on the line between the points (default is False)
        arclen : float (0 < arclen < 1)
            what percentage of a circle the arc between the two points should be (e.g. 0.5 = half circle)

        Returns
        -------
        self

        Examples
        --------

        from manim import *
        import numpy as np

        class Example_FeynManim(Scene):
            def construct(self):

                p1 = np.array([-2,-1,0])
                p2 = np.array([0,2,0])

                p3 = np.array([3,3,0])
                p4 = np.array([0,-2,0])

                example_arc = FeynManim(p2,p1).set_geometry("arc",arclen=np.abs(np.random.random()*0.5)).set_style("photon").build_mobject(YELLOW)
                example_line = FeynManim(p3,p4).set_geometry("line").set_style("gluon").build_mobject(PURPLE)

                self.add(example_arc.mobject)
                self.add(example_line.mobject)

                self.add(Dot(p1),Dot(p2))
                self.add(Dot(p3),Dot(p4))
                
                self.add(Dot(example_arc.center,color=YELLOW))

        """

        if(invert):
            self.vertex1, self.vertex2 = self.vertex2, self.vertex1

        if(geotype=="line"):
            self.geometry = lambda t: self.vertex1 + t*(self.vertex2-self.vertex1)
            self.geolength = np.linalg.norm(self.vertex2-self.vertex1)
            self.center = 0.5*(self.vertex1 + self.vertex2)
            self.momentum_flow_field = lambda pos: (self.vertex2-self.vertex1)/np.linalg.norm(self.vertex2-self.vertex1)

        elif(geotype=="arc"):

            # if(arclen>0.5): # not implemented yet, just do two half circles


            chord = np.linalg.norm(self.vertex2-self.vertex1)
            arclen_angle = arclen*2*np.pi
            self.radius = chord/(2*np.sin(arclen_angle/2))
            height = np.sqrt(self.radius**2 - (chord/2)**2)
            
            self.center = 0.5*(self.vertex1 + self.vertex2) + self._orthonormal_vec((self.vertex2 - self.vertex1)/chord)*height

            # p1p, p2p are coordinates, in which the center of the circle is [0,0]
            p1p = self.vertex1 - self.center
            p2p = self.vertex2 - self.center


            self.geolength = arclen_angle * self.radius

            # shift the start of the circle to vertex1:
            # arccos only valid between 0 and pi 
            # --> unshifted function starts at rightmost point and draws upwards
            # --> if a point lies in the lower half of the circle, take the opposite angle

            if(self.vertex1[1] >= self.center[1]): # starting vertex is in the upper half of the circle
                startingangle = np.arccos(p1p[0]/self.radius) # cosine theorem between unit vector in x direction and p1p
            else:
                startingangle = 2*np.pi - np.arccos(p1p[0]/self.radius)

            self.geometry = lambda t: self.center + self.radius*np.array([ # midpoint of the circle + radius of the circle times circle function
                np.cos(t*arclen_angle + startingangle), 
                np.sin(t*arclen_angle + startingangle), 
                0
            ])

            # multiplied by masking function to be defined above that forces the vector field 0 everywhere but close to the curve
            self.momentum_flow_field = lambda pos: np.array([-pos[1]+self.center[1],pos[0]-self.center[0],0])

        return self



    # # # set styles ----------------------------------------------
    def set_style(self,style,wiggles=False,amp=1):

        """Set the style of the curve

        Parameters
        ----------
        style : str
            either "plain", "fermion", "gluon" or "photon"
        wiggles : int
            number of wiggles/loops on the curve (default is False which forces to match it to the path length)
        amp : float
            width of the gluon and photon styles (default is 1)

        Returns
        -------
        self
        """

        if (style=="plain"):
            self.line_style = lambda t: np.array([0,0,0])

        elif (style=="fermion"):
            self.line_style = lambda t: np.array([0,0,0])
            self.add_arrow()

        elif (style=="gluon"):
            if(wiggles==False):
                wiggles=np.floor(2.5*self.geolength)
            else:
                wiggles += 1

            self.line_style = lambda t: amp*(self.geolength/(2*wiggles)) * ( np.array([0,-1,0]) + np.array([
                np.sin(t*wiggles*2*np.pi), # one wiggle is a full turning of the circle = one loop
                np.cos(t*wiggles*2*np.pi),
                0
            ]))


        elif (style=="photon"):
            if(wiggles==False):
                wiggles=2*np.floor(2.5*self.geolength) # force an even (?) number so start and end of a photon line will always fit together for same winding direction
            else:
                wiggles += 1

            self.line_style = lambda t: amp*(self.geolength/wiggles)*np.array([
                t*0,
                np.sin(t*wiggles*np.pi),
                0
            ])

        return self



    
    # method that builds the flow tube. 
    # - still needs a clever mechanism which scales the density of the lines with the length so that all tubes are per default somewhat equal
    # - still needs to work out cut_start_bleed and cut_end_bleed properly for all cases where geometry is inverted and so on
    def build_flow_tube(self,color=WHITE,delta_x=0.3,delta_y=0.3,stroke_width=2,max_anchors_per_line=30,cut_start_bleed=0,cut_end_bleed=0,invert=False):

        """generates a flow_tube attribute using the current geometry

        TODO: The cut_start_bleed and cut_end_bleed seem to be buggy when the geometry is inverted
        TODO: figure out how to make the flow tube homogeneous between different geometry lengths

        Parameters
        ----------
        cut_start_bleed : float (0 < cut_start_bleed < 1)
            cuts percentage off the start of the flow tube
        cut_end_bleed : float (0 < cut_end_bleed < 1)
            cuts percentage off the end of the flow tube
        invert : bool
            reverses the flow by taking the negative of the momentum vector field (default is False)
        *args : any
            hands over other parameters to Manim's StreamLines method

        Returns
        -------
        self
        """

        if(invert):
            inv=1
        else:
            inv=-1

        c_xmin, c_xmax, c_ymin, c_ymax = self._ext_of_curve(self.geometry)

        self.flow_tube = StreamLines(
            lambda pos: inv * self.momentum_flow_field(pos) * self._masking_function(
                pos, self.geometry, curve_pts=int(10*self.geolength),
                cut_start_bleed=cut_start_bleed, cut_end_bleed=cut_end_bleed
                ), 
            color=color, 
            x_min=c_xmin, x_max=c_xmax, y_min=c_ymin, y_max=c_ymax, 
            delta_x=delta_x, delta_y=delta_y, stroke_width=stroke_width, max_anchors_per_line=max_anchors_per_line)
        return self


    # # # Add Fermion Arrows (and possibly other decorators) -----------------------
    # pos is given in terms of the parameter (t) range
    def add_arrow(self,flip=False,pos=0.5,color=WHITE):

        """generates a Manim Triangle mobject to be used as fermion arrow

        Stores any arrows in a list at self.arrows. Fermion lines add an arrow by default. 

        Parameters
        ----------
        flip : bool
            reverse the arrow direction (default is False)
        pos : float (0 <= pos <= 1)
            places the arrow on the curve in percentage from start to end (default is 0.5)
        color : manim color
            color of the arrow (default is WHITE)

        Returns
        -------
        self
        """

        if(flip):
            flippingangle = np.pi
        else:
            flippingangle = 0

        # angle the object requires to turn the tip towards 0 degree
        triangleangle = np.pi/2

        geometry_gradient = self._derivative(self.geometry)

        # if derivative in x direction + and in y direction +, we are in the first quadrant 
        # if derivative in x direction - and in y direction +, we are in the second quadrant 
        # (starting from 3'clock going counter clockwise)
        if(
            (geometry_gradient(pos)[0] >= 0 and geometry_gradient(pos)[1] >= 0) # x' >= 0 and y' >= 0
            or
            (geometry_gradient(pos)[0] < 0 and geometry_gradient(pos)[1] >= 0) # x' < 0 and y' >= 0
            ):
            rotatingangle = np.arccos(
                np.dot(geometry_gradient(pos),np.array([1,0,0])) / np.linalg.norm(geometry_gradient(pos))
                )
        else:
            rotatingangle = 2*np.pi - np.arccos(
                np.dot(geometry_gradient(pos),np.array([1,0,0])) / np.linalg.norm(geometry_gradient(pos))
                )

        self.arrows.append(Triangle(color=color,fill_opacity=1).scale(0.1).rotate(rotatingangle + flippingangle + triangleangle).move_to(self.geometry(pos)))
        return self

    # sets the self.styled_geometry and the self.mobject element (which are the same unless decorators are present which are included in the mobject)
    def build_mobject(self,color=WHITE):

        """generates a Manim mobject of the curve

        The method first generates self.styled_geometry, which combines self.geometry and self.style. 
        The mobject is then created as a ParametricFunction and accessible at self.mobject
        In case of added arrows, self.mobject is a VGroup of the ParametricFunction and the arrows

        Parameters
        ----------
        color : Manim color
            color of the curve (default is WHITE)

        Returns
        -------
        self
        """


        # gradient of the geometry function
        geometry_gradient = self._derivative(self.geometry)

        # construct a (SO(2)) rotation matrix from the gradient with cos theta = gradx / sqrt(gradx^2 + grady^2); sin theta = grady / sqrt(gradx^2 + grady^2)
        _geometry_gradient_s02 = lambda t: np.array([
            [geometry_gradient(t)[0],-geometry_gradient(t)[1],0],
            [geometry_gradient(t)[1],geometry_gradient(t)[0],0],
            [0,0,0]
            ])/np.linalg.norm(geometry_gradient(t))

        # draw the geometry function plus the gradient rotation matrix times the style function
        self.styled_geometry = lambda t: self.geometry(t) + np.dot(_geometry_gradient_s02(t),self.line_style(t))

        styled_geometry_parametric = ParametricFunction(lambda t: self.styled_geometry(t),t_range=[0,1,0.01],color=color)

        if(len(self.arrows)>0):
            for ar in self.arrows:
                ar.set_color(color)
            self.mobject = VGroup(styled_geometry_parametric,*self.arrows)
        else:
            self.mobject = styled_geometry_parametric

        return self
