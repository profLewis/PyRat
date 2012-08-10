from enthought.tvtk.tools import visual
from enthought.mayavi import mlab


def draw_cyl(coords):
        x_ax = coords['x1'] -coords['x0']
        y_ax = coords['y1'] -coords['y0']
        z_ax = coords['z1'] -coords['z0']
        cyl = visual.Cylinder(pos=(coords['x0'],coords['y0'],coords['z0']),
                          axis=(x_ax,y_ax,z_ax), radius=1, length=10)
        return cyl



        cyl_a_coords = {'y1': 0.0, 'y0': 0.0, 'x0': 0.0, 'x1': 30.0,    'z0': 0.0, 'z1': 0.0}
        cyl_b_coords = {'y1': 0.0, 'y0': 0.0, 'x0': 30.0, 'x1': 630.0, 'z0': 0.0, 'z1': 0.0}
        cyl_c_coords = {'y1': -77.88, 'y0': 0.0, 'x0': 0.0,     'x1': -184.21, 'z0': 0.0, 'z1': 0.0}
        cyl_d_coords = {'y1': 389.41, 'y0': 0.0, 'x0': 0.0,     'x1': -921.06, 'z0': 0.0, 'z1': 0.0}

        cyl_a = {'color': (1.0, 0.0, 0.0), 'coords' :  cyl_a_coords}
        cyl_b = {'color': (1.0, 0.63, 0.04), 'coords' :  cyl_b_coords}
        cyl_c = {'color': (1.0, 0.8, 0.0), 'coords' :  cyl_c_coords}
        cyl_d = {'color': (1.0, 1.0, 0.0), 'coords' :  cyl_d_coords}


        cyls_obj = []
        for cyl in [cyl_a, cyl_b, cyl_c, cyl_d]:
                drawn_cyl = draw_cyl(cyl['coords'])
                drawn_cyl.color = cyl['color']
                cyls_obj.append(drawn_cyl)

