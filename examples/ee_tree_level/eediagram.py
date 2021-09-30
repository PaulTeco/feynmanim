from manim import *
import sys
sys.path.append('../..')
from FeynManim import FeynManim
    

class TreeLevel(Scene):

    def construct(self):

        # generate ee scattering diagram
        topleft = np.array([-3,2,0])
        bottomleft = np.array([-3,-2,0])
        topright = np.array([3,2,0])
        bottomright = np.array([3,-2,0])

        p1 = np.array([-2,0,0])
        vertex1 = Dot(p1)
        p2 = np.array([2,0,0])
        vertex2 = Dot(p2)


        inc1 = FeynManim(bottomleft,p1).set_geometry("line",invert=True).set_style("fermion").build_mobject(WHITE)
        inc1text = Tex(r"$e^-$").next_to(inc1.mobject, LEFT+DOWN)

        inc2 = FeynManim(topleft,p1).set_geometry("line").set_style("fermion").build_mobject(WHITE)
        inc2text = Tex(r"$e^-$").next_to(inc2.mobject, LEFT+UP)

        prop = FeynManim(p1,p2).set_geometry("line").set_style("photon").build_mobject(YELLOW)
        proptext=Tex(r"$\gamma$").scale(1.25).next_to(prop.mobject, UP)

        out1 = FeynManim(p2,bottomright).set_geometry("line",invert=True).set_style("fermion").build_mobject(WHITE)
        out1text = Tex(r"$e^+$").next_to(out1.mobject, RIGHT+DOWN)

        out2 = FeynManim(p2,topright).set_geometry("line").set_style("fermion").build_mobject(WHITE)
        out2text = Tex(r"$e^+$").next_to(out2.mobject, RIGHT+UP)

        self.play(
            LaggedStart(
                Create(inc2text),Create(inc2.mobject),
                Create(inc1.mobject),Create(inc1text),
                Create(vertex1),Create(prop.mobject),Create(proptext),
                Create(out1text),Create(out1.mobject),Create(vertex2),
                Create(out2.mobject),Create(out2text),
                run_time=6,lag_ratio=0.3,rate_func=linear)
        )

        self.wait(2)

        eescattering_elements = VGroup(
            inc1.mobject,inc2.mobject,out1.mobject,out2.mobject,prop.mobject,
            inc1text,inc2text,out1text,out2text,proptext,
            vertex1,vertex2
        )

        eescattering_frame = SurroundingRectangle(eescattering_elements,buff=MED_SMALL_BUFF)

        self.play(Create(eescattering_frame))

        eescattering_symbol = VGroup(eescattering_elements,eescattering_frame)

        eescattering_symbol_target = eescattering_symbol.copy().scale(0.3).shift(5*LEFT+UP)
        
        self.play(ReplacementTransform(eescattering_symbol,eescattering_symbol_target))



        # generate annihilation-creation diagram
        topleft = np.array([-2,3,0])
        bottomleft = np.array([-2,-3,0])
        topright = np.array([2,3,0])
        bottomright = np.array([2,-3,0])

        p1 = np.array([0,-2,0])
        vertex1 = Dot(p1)
        p2 = np.array([0,2,0])
        vertex2 = Dot(p2)

        inc1 = FeynManim(bottomleft,p1).set_geometry("line",invert=True).set_style("fermion").build_mobject(WHITE)
        inc1text = Tex(r"$e^-$").next_to(inc1.mobject, LEFT+DOWN)

        inc2 = FeynManim(p1,bottomright).set_geometry("line",invert=True).set_style("fermion").build_mobject(WHITE)
        inc2text = Tex(r"$e^+$").next_to(inc2.mobject, RIGHT+DOWN)

        prop = FeynManim(p1,p2).set_geometry("line").set_style("photon").build_mobject(YELLOW)
        proptext=Tex(r"$\gamma$").scale(1.25).next_to(prop.mobject, LEFT)

        out1 = FeynManim(topleft,p2).set_geometry("line").set_style("fermion").build_mobject(WHITE)
        out1text = Tex(r"$e^-$").next_to(out1.mobject, LEFT+UP)


        out2 = FeynManim(p2,topright).set_geometry("line").set_style("fermion").build_mobject(WHITE)
        out2text = Tex(r"$e^+$").next_to(out2.mobject, RIGHT+UP)

        self.play(
            LaggedStart(
                Create(inc2text),Create(inc2.mobject),
                Create(inc1.mobject),Create(inc1text),
                Create(vertex1),Create(prop.mobject),Create(proptext),
                Create(out1text),Create(out1.mobject),Create(vertex2),
                Create(out2.mobject),Create(out2text),
                run_time=6,lag_ratio=0.3,rate_func=linear)
        )

        self.wait(2)

        eeannihilation_elements = VGroup(
            inc1.mobject,inc2.mobject,out1.mobject,out2.mobject,prop.mobject,
            inc1text,inc2text,out1text,out2text,proptext,
            vertex1,vertex2
        )

        eeannihilation_frame = SurroundingRectangle(eeannihilation_elements,buff=MED_SMALL_BUFF)

        self.play(Create(eeannihilation_frame))

        eeannihilation_symbol = VGroup(eeannihilation_elements,eeannihilation_frame)

        eeannihilation_symbol_target = eeannihilation_symbol.copy().scale(0.3).shift(5*RIGHT+UP)
        
        self.play(ReplacementTransform(eeannihilation_symbol,eeannihilation_symbol_target))

        self.wait(0.5)

        plus = Tex(r"$+$").scale(2).shift(UP)
        eescattering_symbol_target2 = eescattering_symbol_target.copy().shift(2*RIGHT).scale(1.4)
        eeannihilation_symbol_target2 = eeannihilation_symbol_target.copy().shift(2*LEFT).scale(1.4)


        self.play(
            ReplacementTransform(eescattering_symbol_target,eescattering_symbol_target2),
            Create(plus),
            ReplacementTransform(eeannihilation_symbol_target,eeannihilation_symbol_target2)
        )

        tree_level_caption = Text("Tree Level",color=YELLOW).scale(2).shift(2*DOWN)

        self.play(FadeIn(tree_level_caption))

        self.wait(2)


