from manim import *
import sys
sys.path.append('../..')
from FeynManim import FeynManim
    


# instantiate diagram with defined "diagram reading order"
gluonpenguin_full = {
    "inctop": {"order":0},
    "vertex1": {"order":1},
    "loopdownleft": {"order":2},
    "vertex3": {"order":3},
    "loopdownright": {"order":4},
    "vertex2": {"order":5},
    "outtop": {"order":6},
    "outbottom": {"order":7},
    "vertex4": {"order":8},
    "incbottom": {"order":9},
    "wprop": {"order":10},
    "gluonprop": {"order":11},
}


# add feynman diagram elements

# coordinate definitions first
NW_corner = [-3,3,0]
NE_corner = [3,3,0]
SE_corner = [3,-2,0]
SW_corner = [-3,-2,0]

gluonpenguin_full["vertex1"]["manim"] = Dot([-2,2,0])
gluonpenguin_full["vertex2"]["manim"] = Dot([2,2,0])
gluonpenguin_full["vertex3"]["manim"] = Dot([0,0,0])
gluonpenguin_full["vertex4"]["manim"] = Dot([0,-2,0])

# # then the rest in "diagram reading order"
gluonpenguin_full["inctop"]["feyn"] = FeynManim(NW_corner,gluonpenguin_full["vertex1"]["manim"].get_center()).set_geometry("line").set_style("fermion").build_mobject()
gluonpenguin_full["loopdownleft"]["feyn"] = FeynManim(gluonpenguin_full["vertex1"]["manim"].get_center(),gluonpenguin_full["vertex3"]["manim"].get_center()).set_geometry("arc",arclen=0.25).set_style("fermion").build_mobject()
gluonpenguin_full["loopdownright"]["feyn"] = FeynManim(gluonpenguin_full["vertex3"]["manim"].get_center(),gluonpenguin_full["vertex2"]["manim"].get_center()).set_geometry("arc",arclen=0.25).set_style("fermion").build_mobject()
gluonpenguin_full["outtop"]["feyn"] = FeynManim(gluonpenguin_full["vertex2"]["manim"].get_center(),NE_corner).set_geometry("line").set_style("fermion").build_mobject()
gluonpenguin_full["outbottom"]["feyn"] = FeynManim(gluonpenguin_full["vertex4"]["manim"].get_center(),SE_corner).set_geometry("line").set_style("fermion").build_mobject()
gluonpenguin_full["incbottom"]["feyn"] = FeynManim(SW_corner,gluonpenguin_full["vertex4"]["manim"].get_center()).set_geometry("line").set_style("fermion").build_mobject()
gluonpenguin_full["wprop"]["feyn"] = FeynManim(gluonpenguin_full["vertex1"]["manim"].get_center(),gluonpenguin_full["vertex2"]["manim"].get_center()).set_geometry("arc",arclen=0.125,invert=True).set_style("photon").build_mobject(color=YELLOW)
gluonpenguin_full["gluonprop"]["feyn"] = FeynManim(gluonpenguin_full["vertex3"]["manim"].get_center(),gluonpenguin_full["vertex4"]["manim"].get_center()).set_geometry("line").set_style("gluon").build_mobject(color=YELLOW)

# # set manim objects for quicker reference
gluonpenguin_full["inctop"]["manim"] = gluonpenguin_full["inctop"]["feyn"].mobject
gluonpenguin_full["loopdownleft"]["manim"] = gluonpenguin_full["loopdownleft"]["feyn"].mobject
gluonpenguin_full["loopdownright"]["manim"] = gluonpenguin_full["loopdownright"]["feyn"].mobject
gluonpenguin_full["outtop"]["manim"] = gluonpenguin_full["outtop"]["feyn"].mobject
gluonpenguin_full["outbottom"]["manim"] = gluonpenguin_full["outbottom"]["feyn"].mobject
gluonpenguin_full["incbottom"]["manim"] = gluonpenguin_full["incbottom"]["feyn"].mobject
gluonpenguin_full["wprop"]["manim"] = gluonpenguin_full["wprop"]["feyn"].mobject
gluonpenguin_full["gluonprop"]["manim"] = gluonpenguin_full["gluonprop"]["feyn"].mobject




# add labels to the diagram elements

# in "diagram reading order"
gluonpenguin_full["inctop"]["label"] = Tex(r"$\bar{s}$").next_to(gluonpenguin_full["inctop"]["manim"], LEFT)
gluonpenguin_full["vertex1"]["label"] = False
gluonpenguin_full["loopdownleft"]["label"] = Tex(r"$\bar{q} = \bar{u},\bar{c},\bar{t}$").next_to(gluonpenguin_full["vertex3"]["manim"].get_center(),2*UP)
gluonpenguin_full["vertex3"]["label"] = False
gluonpenguin_full["loopdownright"]["label"] = False
gluonpenguin_full["vertex2"]["label"] = False
gluonpenguin_full["outtop"]["label"] = Tex(r"$\bar{d}$").next_to(gluonpenguin_full["outtop"]["manim"], RIGHT)
gluonpenguin_full["outbottom"]["label"] = Tex(r"$u$").next_to(gluonpenguin_full["outbottom"]["manim"], RIGHT)
gluonpenguin_full["vertex4"]["label"] = False
gluonpenguin_full["incbottom"]["label"] = Tex(r"$u$").next_to(gluonpenguin_full["incbottom"]["manim"], LEFT)
gluonpenguin_full["wprop"]["label"]=Tex(r"$W$").next_to(gluonpenguin_full["wprop"]["manim"], UP)
gluonpenguin_full["gluonprop"]["label"] = Tex(r"$g$").next_to(gluonpenguin_full["gluonprop"]["manim"],RIGHT)





# add equation terms which the diagram elements represent

gluonpenguin_full["inctop"]["eq"] = r"\bar{s}_i"
gluonpenguin_full["vertex1"]["eq"] = r"\left(-i \eta \frac{g}{\sqrt{2}} \gamma^{\mu} P_L V_{us}^{*} \right)"
gluonpenguin_full["loopdownleft"]["eq"] = r"\left( i\frac{\left( \slashed{p}_1 - m_q \right)}{p_1^2 + m_q^2 - i \epsilon} \right)"
gluonpenguin_full["vertex3"]["eq"] = r"\left(-i \eta_s g_s \gamma^{\rho} T^a_{ij} \right)"
gluonpenguin_full["loopdownright"]["eq"] = r"\left( i\frac{\left( \slashed{p}_2 - m_q \right)}{p_2^2 + m_q^2 - i \epsilon} \right)\\"
gluonpenguin_full["vertex2"]["eq"] = r"&\left(-i \eta \frac{g}{\sqrt{2}} \gamma^{\nu} P_L V_{ud} \right)"
gluonpenguin_full["outtop"]["eq"] = r"d_j"
gluonpenguin_full["outbottom"]["eq"] = r"\bar{u_k}"
gluonpenguin_full["vertex4"]["eq"] = r"\left(-i \eta_s g_s \gamma^{\sigma} T^b_{lk} \right)"
gluonpenguin_full["incbottom"]["eq"] = r"u_l\\"
gluonpenguin_full["wprop"]["eq"] = r"&\left(-i \frac{1}{p^2_3 - m_W^2 + i \epsilon} \left[ g_{\mu \nu} - \left( 1 - \xi_W \right) \frac{p_{3\mu}p_{3\nu}}{p_3^2 - \xi_W m_W^2}\right]\right)"
gluonpenguin_full["gluonprop"]["eq"] = r"\left(-i \delta_{ab} \frac{1}{p_4^2 + i \epsilon} \left[ g_{\rho \sigma} - \left( 1 - \xi_G \right) \frac{p_{4\rho}p_{4\sigma}}{p_4^2}\right]\right)"


# lists of the respective parts for easier drawing
feyn_full = [gluonpenguin_full[partname]["manim"] for partname in gluonpenguin_full.keys() if gluonpenguin_full[partname]["manim"]]
label_full = [gluonpenguin_full[partname]["label"] for partname in gluonpenguin_full.keys() if gluonpenguin_full[partname]["label"]]
eq_full = [gluonpenguin_full[partname]["eq"] for partname in gluonpenguin_full.keys() if gluonpenguin_full[partname]["eq"]]

# list that groups diagram parts and their labels (if one exists) together
feyn_label_full = [
    VGroup(*[el for el in [gluonpenguin_full[partname]["manim"], gluonpenguin_full[partname]["label"]] if el]) 
    for partname in gluonpenguin_full.keys()
    if gluonpenguin_full[partname]["manim"]
    ]

# list that sorts create objects from left to right as lists. loop over the list elements and use self.play(*element)
feyn_label_create_timeorder = [[Create(part["manim"]),Create(part["label"])] if part["label"] else [Create(part["manim"])]
for part in sorted([gluonpenguin_full[el] for el in gluonpenguin_full.keys()], key=lambda k: k["manim"].get_center()[0])
if part["manim"]
]



# set up the equations

slashed = TexTemplate()
slashed.add_to_preamble(r"\usepackage{slashed}")

iM_eq = r"i \mathcal{M}=&"

full_eq_tex=MathTex(
    iM_eq, *eq_full, tex_template=slashed
).shift(2*DOWN).scale(0.6)

full_eq_tex.set_color_by_tex(gluonpenguin_full["wprop"]["eq"], YELLOW)
full_eq_tex.set_color_by_tex(gluonpenguin_full["gluonprop"]["eq"], YELLOW)







##############################################################################################################
#
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # SCENES # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#
#
##############################################################################################################









# Test Scene to try various stuff
class Test(Scene):

    def construct(self):
        p1 = np.array([1,1,0])
        p2 = np.array([0,2,0])
        l1 = FeynManim(p1,p2).set_geometry("arc").set_style("gluon").build_mobject(YELLOW)

        self.add(l1.mobject,Dot(p1),Dot(p2))



# Drawing the diagram and writing down the Feynman Rules (c.f. https://arxiv.org/abs/1209.6213)
class GluonPenguin(Scene):

    def construct(self):

        # # # Different Steps in Pictures

        # # Step 1 Drawing the Diagram with labels
        # # self.add(VGroup(*feyn_full))
        # # self.add(VGroup(*label_full))


        # # Step 2 Shrinking the Diagram and writing the full equation
        # self.add(VGroup(*feyn_full).scale(0.6).shift(UP))
        # self.add(VGroup(*label_full).scale(0.6).shift(UP))
        # self.add(full_eq_tex)

        # # Step 3: Removing the Equation and rescaling the Diagram without labels
        # # self.add(VGroup(*feyn_full))

        # # # --------------------------------------------

        # # # Animation

        # Draw the Diagram
        for part in feyn_label_create_timeorder:
            self.play(*part)

        # Shrink the diagram and the labels to make space for equations
        feyn_full_shrink = VGroup(*feyn_full).copy().scale(0.6).shift(UP)
        label_full_shrink = VGroup(*label_full).copy().scale(0.6).shift(UP)
        self.play(ReplacementTransform(VGroup(*feyn_full),feyn_full_shrink),ReplacementTransform(VGroup(*label_full),label_full_shrink))
        self.wait()

        # Define frameboxes
        frameboxes_start = []
        frameboxes_finish = []
        for i in range(len(feyn_label_full)):
            frameboxes_start += [SurroundingRectangle(feyn_label_full[i], buff = .2, color=BLUE)]
            frameboxes_finish += [SurroundingRectangle(full_eq_tex[i+1], buff = .1, color=BLUE)]

        self.play(Create(full_eq_tex[0]))

        # animate frameboxes and build up the equation
        for i in range(len(feyn_label_full)):
            self.play(Create(frameboxes_start[i]))
            self.wait()
            self.play(ReplacementTransform(frameboxes_start[i],frameboxes_finish[i]))
            self.play(Write(full_eq_tex[i+1]))
            self.wait()
            self.play(FadeOut(frameboxes_finish[i]))
            self.wait()

        # enlarge again to original size without equation
        self.play(FadeOut(label_full_shrink),FadeOut(full_eq_tex))
        self.play(ReplacementTransform(feyn_full_shrink,feyn_full_shrink.copy().shift(DOWN).scale(4/3)))

        self.wait(2)


# Momentum Flow Animation
class MomentumFlow(Scene):

    def construct(self):
        
        diagram_group = VGroup(*feyn_full)


        self.add(diagram_group)
        
        diagram_group_fade = diagram_group.copy().set_stroke(opacity=0.2)

        self.play(ReplacementTransform(diagram_group,diagram_group_fade))

        gluonpenguin_full["inctop"]["feyn"].build_flow_tube(color=YELLOW,cut_end_bleed=0.3,invert=True)
        gluonpenguin_full["loopdownleft"]["feyn"].build_flow_tube(color=YELLOW,cut_end_bleed=0.2,invert=True)
        gluonpenguin_full["loopdownright"]["feyn"].build_flow_tube(color=YELLOW,cut_end_bleed=0.2,invert=True)
        gluonpenguin_full["outtop"]["feyn"].build_flow_tube(color=YELLOW,invert=True)

        gluonpenguin_full["outbottom"]["feyn"].build_flow_tube(color=YELLOW,invert=True)
        gluonpenguin_full["incbottom"]["feyn"].build_flow_tube(color=YELLOW,invert=True)

        gluonpenguin_full["gluonprop"]["feyn"].build_flow_tube(color=YELLOW,invert=True)
        gluonpenguin_full["wprop"]["feyn"].build_flow_tube(color=YELLOW,invert=True)

        # Animate Momentum Flow
        stream_lines = [
            gluonpenguin_full["inctop"]["feyn"].flow_tube,
            gluonpenguin_full["loopdownleft"]["feyn"].flow_tube,
            gluonpenguin_full["loopdownright"]["feyn"].flow_tube,
            gluonpenguin_full["outtop"]["feyn"].flow_tube,
            gluonpenguin_full["outbottom"]["feyn"].flow_tube,
            gluonpenguin_full["incbottom"]["feyn"].flow_tube,
            gluonpenguin_full["gluonprop"]["feyn"].flow_tube,
            gluonpenguin_full["wprop"]["feyn"].flow_tube
            ]

        self.add(*stream_lines)

        flow_animation_time = 1

        for s in stream_lines:
            s.start_animation(warm_up=True, flow_speed=1.5)
        
        self.wait(flow_animation_time * stream_lines[0].virtual_time / stream_lines[0].flow_speed)


# Generate Momentum Flow Picture
class MomentumArrows(Scene):

    def construct(self):
        
        diagram_group = VGroup(*feyn_full)


        self.add(diagram_group)

        arrow_inctop = Arrow(
            start=(LEFT+UP)/(2*np.sqrt(2)), end=(RIGHT+DOWN)/(2*np.sqrt(2)), buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["inctop"]["manim"], LEFT)


        arrow_loopdownleft = Arrow(
            start=(LEFT+UP)/(2*np.sqrt(2)), end=(RIGHT+DOWN)/(2*np.sqrt(2)), buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["loopdownleft"]["manim"], LEFT+DOWN).shift(3*(RIGHT+UP)/4)

        arrow_loopdownright = Arrow(
            start=(LEFT+DOWN)/(2*np.sqrt(2)), end=(RIGHT+UP)/(2*np.sqrt(2)), buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["loopdownright"]["manim"], (RIGHT+DOWN)).shift(3*(LEFT+UP)/4)

        arrow_outtop = Arrow(
            start=(LEFT+DOWN)/(2*np.sqrt(2)), end=(RIGHT+UP)/(2*np.sqrt(2)), buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["outtop"]["manim"], RIGHT)
        
        arrow_wprop = Arrow(
            start=(RIGHT)/2, end=(LEFT)/2, buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["wprop"]["manim"], UP)

        arrow_gluonprop = Arrow(
            start=(UP)/2, end=(DOWN)/2, buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["gluonprop"]["manim"], LEFT)


        arrow_outbottom = Arrow(
            start=RIGHT/2, end=LEFT/2, buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["outbottom"]["manim"], DOWN)

        arrow_incbottom = Arrow(
            start=RIGHT/2, end=LEFT/2, buff=0,
            max_tip_length_to_length_ratio=0.15,max_stroke_width_to_length_ratio=3
            ).next_to(gluonpenguin_full["incbottom"]["manim"], DOWN)


        self.add(arrow_inctop,arrow_loopdownleft,arrow_loopdownright,arrow_outtop,arrow_wprop,arrow_gluonprop,arrow_outbottom,arrow_incbottom)










