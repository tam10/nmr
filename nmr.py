# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 11:50:10 2017

@author: Tristan Mackenzie

    QNMR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    QNMR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with QNMR.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import scipy.ndimage as ndi
import scipy.interpolate as interpolate

import sys
if sys.version_info[0] == 3:
    import tkinter as tk
    import tkinter.messagebox as msgbox
    from tkinter.filedialog import askopenfilename, asksaveasfilename
else:
    import Tkinter as tk
    import tkMessageBox as msgbox
    from tkFileDialog import askopenfilename, asksaveasfilename

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib import pyplot as plt

class GUI():
    
    def __init__(self):
        
        self.root = tk.Tk()
        self.root.title("NMR GUI")
        self.root.resizable(0,0)
        
        sunken = dict(height = 2, bd = 1, relief = "sunken")        
        
        self.figure = plt.figure(figsize = (10, 5))
        self.ax     = self.figure.add_subplot(111)
        self.ax.invert_xaxis()
        
        self.peak       = None
        self.peaks      = []
        self.splitting  = None
        self.splittings = []
        
        fs  = self.frames            = {}
        cs  = self.canvases          = {}
        ls  = self.labels            = {}
        mes = self.machine_entries   = {}
        pes = self.peak_entries      = {}
        ses = self.splitting_entries = {}
        bs  = self.buttons           = {}
        ms  = self.optionmenus       = {}
        
        fs ["machine"]               = _add_frame(dict(master=self.root, text="Machine", **sunken), gk('000055news'))
        ls ["frequency"]             = _add_label(fs["machine"], {"text": "Operating Frequency (MHz):"}, gk('00w'))
        mes["machine_frequency_mhz"] = _add_entry(fs["machine"], "", {}, gk('010300'))
        ls ["noise"]                 = _add_label(fs["machine"], {"text": "Noise:"}, gk('10w'))
        mes["noise"]                 = _add_entry(fs["machine"], "", {}, gk('110300'))
        ls ["resolution"]            = _add_label(fs["machine"], {"text": "Resolution (ppm):"}, gk('20w'))
        mes["resolution"]            = _add_entry(fs["machine"], "", {}, gk('210300'))
        ls ["min_x"]                 = _add_label(fs["machine"], {"text": "Range (ppm):"}, gk('30w'))
        mes["min_x"]                 = _add_entry(fs["machine"], "", {}, gk('31w'), {"width": 3})
        ls ["max_x"]                 = _add_label(fs["machine"], {"text": "to:"}, gk('32w'))
        mes["max_x"]                 = _add_entry(fs["machine"], "", {}, gk('33w'), {"width": 3})
        
        fs["peaks"]            = _add_frame(dict(master=self.root, text="Peaks", **sunken), gk('100055news'))
        bs["add_peak"]         = _add_button(fs["peaks"], {"text": "Add Peak"}, gk('000055w'), {"<Button-1>": self.add_peak})
        bs["remove_peak"]      = _add_button(fs["peaks"], {"text": "Remove Peak"}, gk('010055w'), {"<Button-1>": self.remove_peak})
        ls["peaks"]            = _add_label(fs["peaks"], {"text": "Peaks:"}, gk('10w'))
        ms["peaks"], self.peak_string = _add_optionmenu(fs["peaks"], " ", [" "], {"command": self._update_peak_om}, gk('1103ew'), {"width": 10})
        
        ls ["shift"]    = _add_label(fs["peaks"], {"text": "Shift (ppm):"}, gk('20w'))
        pes["shift"]    = _add_entry(fs["peaks"], "", {}, gk('21w'), attach_func=self._set_peak_string)
        ls ["p_nuclei"] = _add_label(fs["peaks"], {"text": "Nuclei:"}, gk('22w'))
        pes["nuclei"]   = _add_entry(fs["peaks"], "", {}, gk('23w'), attach_func=self._set_peak_string)
        ls ["hwhm"]     = _add_label(fs["peaks"], {"text": "Half Width Half Maximum (ppm):"}, gk('3003w'))
        pes["hwhm"]     = _add_entry(fs["peaks"], "", {}, gk('33w'), attach_func=self._set_peak_string)
        
        fs["splittings"]       = _add_frame(dict(master=fs["peaks"], text="Splitting Nuclei", **sunken), gk('400455news'))
        bs["add_splitting"]    = _add_button(fs["splittings"], {"text": "Add Splitting"}, gk('000055w'), {"<Button-1>": self.add_splitting})
        bs["remove_splitting"] = _add_button(fs["splittings"], {"text": "Remove Splitting"}, gk('010055w'), {"<Button-1>": self.remove_splitting})
        ls["splittings"]       = _add_label(fs["splittings"], {"text": "Splittings:"}, gk('10w'))
        ms["splittings"], self.splitting_string = _add_optionmenu(fs["splittings"], " ", [" "], {}, gk('1103ew'), {"width": 10})
        
        ls ["coupling"]     = _add_label(fs["splittings"], {"text": "J Coupling (Hz):"}, gk('20w'))
        ses["coupling"]     = _add_entry(fs["splittings"], "", {}, gk('21w'), attach_func=self._set_splitting_string)
        ls ["s_nuclei"]     = _add_label(fs["splittings"], {"text": "Nuclei:"}, gk('22w'))
        ses["nuclei"]       = _add_entry(fs["splittings"], "", {}, gk('23w'), attach_func=self._set_splitting_string)
        ls ["spin"]         = _add_label(fs["splittings"], {"text": "Spin:"}, gk('30w'))
        ses["spin"]         = _add_entry(fs["splittings"], "", {}, gk('31w'), attach_func=self._set_splitting_string)
        ls ["abundance"]    = _add_label(fs["splittings"], {"text": "Abundance:"}, gk('32w'))
        ses["abundance"]    = _add_entry(fs["splittings"], "", {}, gk('33w'), attach_func=self._set_splitting_string)
        
        fs["controls"] = _add_frame(dict(master=self.root, text="Controls", **sunken), gk('200055news'))
        bs["update"]   = _add_button(fs["controls"], {"text": "Update"}, gk('000055ew') ,{"<Button-1>": self.update})
        bs["parse"]    = _add_button(fs["controls"], {"text": "From .log"}, gk('010055ew') ,{"<Button-1>": self.parse})
        bs["export"]   = _add_button(fs["controls"], {"text": "Export Data"}, gk('020055ew') ,{"<Button-1>": self.export})
        
        fs["plot"]  = _add_frame(dict(master=self.root, text="Plot", **sunken), gk('012055news'))
        cs["plot"]  = _add_mpl_canvas(fs["plot"], self.figure, gk('00'))
        cs["plot"].get_tk_widget().grid(row=0, column=0)
        
        fs["toolbar"] = _add_frame(dict(master=self.root, text="", **sunken), gk('210055news'))
        
        self.toolbar = NavigationToolbar2TkAgg(cs["plot"], fs["toolbar"])
        self.toolbar.grid(row=0, column=0)
        
        self._add_nmr()
        self._add_nmr_parser()
            
        self.root.protocol("WM_DELETE_WINDOW", self._cl)
        self.root.mainloop()
        
    def _add_nmr(self, *args):
        
        self.nmr = NMR()
        
        for key, entry in self.machine_entries.items():
            entry.configure(textvariable=getattr(self.nmr, key), state=tk.NORMAL)
        
    def _add_nmr_parser(self, *args):
        self.nmr_parser = NMRParser()
        self.nmr_parser._ask_spin_abundance = self._ask_spin_abundance
        
    def _cl(self):            
        plt.close('all')
        self.root.destroy()
        
    def _set_peak(self, peak, *args):
        self.peak = peak
        self._peak_changed()
        
        try:
            splitting = self.peak.splittings[0]
        except (IndexError, AttributeError):
            splitting = None
            
        self.splittings = self.peak.splittings
        self._set_splitting(splitting)

    def _set_peak_string(self, *args):
        self.peak_string.set(repr(self.peak))   
        self._update_peak_om()     
        
    def _peak_changed(self, *args):
        self._set_peak_string()
        self._update_peak_entries()
        
    def _update_peak_om(self):
        
        om = self.optionmenus["peaks"]
        menu = om['menu']
        menu.delete(0, tk.END)
        
        for peak in self.peaks:
            string = repr(peak)
            menu.add_command(label = string, command = lambda value=peak: self._set_peak(value))
            
    def add_peak(self, *args):
        
        peak = Peak(self.nmr, 1, 7)
            
        self.peaks.append(peak)
        self.nmr.peaks.append(peak)
        self._set_peak(peak)
        
    def remove_peak(self, *args):
        
        peak = self.peak
        
        self.nmr.peaks.remove(peak)
        self.peaks.remove(peak)
        
        try:
            peak = self.peaks[0]
        except IndexError:
            peak = None
        
        self._set_peak(peak)
            
    def _update_peak_entries(self, *args):
        peak = self.peak
        if peak:
            for key, entry in self.peak_entries.items():
                entry.configure(textvariable=getattr(peak, key), state=tk.NORMAL)
        else:
            for key, entry in self.peak_entries.items():
                entry.configure(textvariable=tk.StringVar(value=""), state=tk.DISABLED)
            
    def _set_splitting(self, splitting, *args):
        self.splitting = splitting
        self._splitting_changed()

    def _set_splitting_string(self, *args):
        self.splitting_string.set(repr(self.splitting))   
        self._update_splitting_om()     
        
    def _splitting_changed(self, *args):
        self._set_splitting_string()
        self._update_splitting_entries()
        
    def _update_splitting_om(self):
        
        om = self.optionmenus["splittings"]
        menu = om['menu']
        menu.delete(0, tk.END)
        
        for splitting in self.splittings:
            string = repr(splitting)
            menu.add_command(label = string, command = lambda value=splitting: self._set_splitting(value))
            
    def add_splitting(self, *args):
        
        splitting = Splitting(0.5, 1, 20, 1)
            
        self.splittings.append(splitting)
        self.peak.splittings.append(splitting)
        self._set_splitting(splitting)
        
    def remove_splitting(self, *args):
        
        s0 = self.splitting
        
        for i, s1 in enumerate(self.peak.splittings):
            if s0 == s1:
                del self.peak.splittings[i]
                break
            
        for i, s1 in enumerate(self.splittings):
            if s0 == s1:
                del self.peak.splittings[i]
                break
        
        try:
            splitting = self.splittings[0]
        except IndexError:
            splitting = None
        
        self._set_splitting(splitting)
            
    def _update_splitting_entries(self, *args):
        splitting = self.splitting
        if splitting:
            for key, entry in self.splitting_entries.items():
                entry.configure(textvariable=getattr(splitting, key), state=tk.NORMAL)
        else:
            for key, entry in self.splitting_entries.items():
                entry.configure(textvariable=tk.StringVar(value=""), state=tk.DISABLED)
            
    def _ask_spin_abundance(self, element):
        
        while True:
            sp = EntryPopup(self, "Input nuclear spin for element {}:".format(element))
            sp.root.wait_window()
            spin = sp.value
            try:
                spin = float(spin)
                if spin % 0.5 == 0 and spin >= 0:
                    break
            except:
                pass
            
            msgbox.showerror("Error", "Spin must be positive half-integer")
        
        while True:
            sp = EntryPopup(self, "Input abundance (0-1) for element {}:".format(element))
            sp.root.wait_window()
            abundance = sp.value
            try:
                abundance = float(abundance)
                if 0 < abundance < 1:
                    break
            except:
                pass
            
            msgbox.showerror("Error", "Abundance must be between 0 and 1")
        
    def update(self, *args):    
        
        xs, ys = self.nmr.get_plot()
        
        min_x  = float(self.nmr.min_x.get())
        max_x  = float(self.nmr.max_x.get())
        
        self.ax.clear()
        self.ax.plot(xs, ys)
        self.ax.set_xlim(min_x, max_x)
        self.ax.set_xlabel("Chemical Shift (ppm)")
        self.ax.yaxis.set_visible(False)
        self.ax.invert_xaxis()
        self.figure.tight_layout()
        self.canvases["plot"].draw()
            
    def parse(self, *args):
        
        fn = askopenfilename(filetypes = (("Gaussian Log File", "*.log"), ("All Files", "*.*")))
        
        self.nmr_parser.parse(fn)
        
        gp = LoadGaussianPopup(self, self.nmr_parser)
        gp.root.wait_window()
        shifts = []
        
        try:
            gaussian_nmr_list = self.gaussian_nmr_list
            self.nmr.peaks = []
            self.peaks = []
            self.splittings = []
            for shift, splittings in gaussian_nmr_list:
                
                shifts.append(shift)
                peak = Peak(self.nmr, 1, shift)
                
                self.peaks.append(peak)
                self.nmr.peaks.append(peak)
                    
                for coupling, spin, degeneracy in splittings:
                    
                    splitting = Splitting(spin, degeneracy, coupling, 1)
                    peak.splittings.append(splitting)
                    
            for splitting in peak.splittings:
                self.splittings.append(splitting)
                    
            self._set_peak(peak)
            self._set_splitting(splitting)
            
            
            min_x = min(shifts)
            max_x = max(shifts)
            diff = max_x - min_x
            
            self.nmr.min_x.set(round(min_x - 0.2 * diff) - 1)
            self.nmr.max_x.set(round(max_x + 0.2 * diff) + 1)
            
            self.update()
            
        except:
            
            msgbox.showerror("Error", "Could not load Gaussian .log File")
            raise 

    def export(self, *args):
        
        try:
            line = self.ax.lines[0]
        except IndexError:
            msgbox.showerror("No Data", "No data to export!")
            return
            
        data = line.get_xydata()
        
        fn = asksaveasfilename(filetypes = [("CSV Files", "*.csv")])
        
        with open(fn, "w") as f:
            for row in data:
                f.write("{},{}\n".format(*row))

class NMR():

    def __init__(self):
        self.machine_frequency_mhz = tk.StringVar(value='400')
        self.peaks = []
        self.resolution = tk.StringVar(value='0.01')
        self._epsilon   = tk.StringVar(value='1e-6')
        self.noise      = tk.StringVar(value='0.1')
        self.min_x      = tk.StringVar(value='0')
        self.max_x      = tk.StringVar(value='12')
        
    def get_plot(self):
        
        min_x = float(self.min_x.get())
        max_x = float(self.max_x.get())
        res   = float(self.resolution.get())
        noise = float(self.noise.get())
            
        xs = np.arange(min_x, max_x + res, res)
        ys = np.random.random(len(xs)) * noise
        
        for i, peak in enumerate(self.peaks):
            
            p_xs, p_ys = peak.generate(min_x, max_x)
            
            p_y_ints = interpolate.griddata(p_xs, p_ys, xs, method='linear')
            
            ys += p_y_ints
            
        return xs, ys
        
    def __repr__(self):
        
        return "NMR(freq={}, resolution={}, noise={}, min_x={}, max_x={}".format(
            self.machine_frequency_mhz.get(),
            self.resolution.get(),
            self.noise.get(),
            self.min_x.get(),
            self.max_x.get()
        )
    
class Peak():
    
    def __init__(self, parent, nuclei, shift, hwhm=0.01):
        
        self.nuclei = tk.StringVar(value=nuclei)
        self.shift  = tk.StringVar(value=shift)
        self.hwhm   = tk.StringVar(value=hwhm)
        self.parent = parent
        self.splittings = []
        
    def cauchy(self, min_x, max_x):
        
        res    = float(self.parent.resolution.get())
        hwhm   = float(self.hwhm.get())
        mf     = float(self.parent.machine_frequency_mhz.get())
        nuclei = float(self.nuclei.get())
        shift  = float(self.shift.get())
        
        #Extend x domain to include off-chart contributions to splitting + FWHM
        max_split = 2 * hwhm
        for S in self.splittings:
            nuclei = int(S.nuclei.get())
            spin = float(S.spin.get())
            coupling = float(S.coupling.get())
            
            max_split += (coupling * (nuclei * spin + 1) / mf)
            
        min_x -= (round(max_split / res)) * res
        max_x += (round(max_split / res)) * res
        
        xs = np.arange(min_x, max_x + res, res)
        ys = []
        
        for x in xs:
            ys.append((nuclei / (np.pi * hwhm * (1 + ((x - shift) / hwhm) ** 2))))
        
        return xs, ys
        
    def generate(self, min_x, max_x):
        
        res     = float(self.parent.resolution.get())
        epsilon = float(self.parent._epsilon.get())
        mf      = float(self.parent.machine_frequency_mhz.get())
            
        xs, ys = self.cauchy(min_x, max_x)
        
        if len(xs) == 0:
            return [], []
            
        for S in self.splittings:
            
            nuclei = int(S.nuclei.get())
            spin = float(S.spin.get())
            coupling = float(S.coupling.get())
            
            s = list(S.get_splitting())
            
            j_split = float(coupling) / mf
            
            max_j = (nuclei * spin) * j_split
            
            conv_xs = np.arange(- max_j, max_j + res, res)
            conv_ys = []
            
            j = - max_j
            
            for i, conv_x in enumerate(conv_xs):
                
                if j - conv_x <= epsilon:
                    conv_ys.append(s.pop(0))
                    j += j_split * 0.5
                else:
                    conv_ys.append(0.0)

                
            ys = ndi.convolve1d(ys, conv_ys)
            
            
        return xs, np.array(ys)
             
    def __repr__(self):
        
        return "Peak(nuclei={}, shift={:.3f}, hwhm={:.3f})".format(int(self.nuclei.get()), float(self.shift.get()), float(self.hwhm.get()))
        
class Splitting():
    
    def __init__(self, spin, nuclei, coupling, abundance):
        
        self.spin      = tk.StringVar(value=spin)
        self.nuclei    = tk.StringVar(value=nuclei)
        self.coupling  = tk.StringVar(value=coupling)
        self.abundance = tk.StringVar(value=abundance)
        
    def get_splitting(self):
        
        abundance = float(self.abundance.get())
        
        row = self.pascal()
        norm = sum(row)
        
        row *= abundance / norm
        
        mid = (len(row) - 1) / 2
        
        row[mid] += 1 - abundance
        
        return row
        
    def pascal(self):

        spin   = float(self.spin.get())
        nuclei = int(self.nuclei.get())
        
        if not spin % 0.5 == 0:
            raise ValueError("Spin must be divisible by 0.5")
        
        #Number of elements
        n = int(4 * spin * nuclei + 1)

        prev_row = [1 if i == 2 * spin * nuclei else 0 for i in range(n)]

        for nucleus in range(nuclei):
            row = []

            for i, element in enumerate(range(n)):
                v = 0
                for p_i, p_element in enumerate(prev_row):
                    if abs(p_i - i) <= 2 * spin and (p_i - i) % 2 == 2 * spin % 2:
                        v += p_element

                row.append(float(v))

            prev_row = row

        return np.array(row)
        
    def __repr__(self):
        
        return "Splitting(spin={:.1f}, nuclei={}, coupling={:.3f}, abundance={:.3%})".format(
            float(self.spin.get()),
            int(self.nuclei.get()),
            float(self.coupling.get()),
            float(self.abundance.get())
        )
        
    def __eq__(self, other):
        if isinstance(other, Splitting):
            for a in ["spin", "nuclei", "coupling", "abundance"]:
                if getattr(self, a).get() != getattr(other, a).get():
                    return False
            return True
        else:
            return False

class NMRParser():
    def __init__(self):
        
        self.peak_dict = {}
        self.coupling_degeneracy_threshold = 1
        self.spin_dict = {
            "H" : [0.5, 1],
            "C" : [0.5, 0.011],
            "N" : [0.5, 0.00365],
            "O" : [0, 0],
            "S" : [1.5, 0.0076],
            "Si": [0.5, 0.047]
        }
        
    def parse(self, fn):
        with open(fn, "r") as f:
            lines = f.readlines()

        ln = 0
        n_ln = len(lines)
        n_ats = 0
        elements = []
        shifts = []
        spins = []


        while ln < n_ln:
            line = lines[ln]
            if n_ats == 0:
                if line.strip() in ["Input orientation:", "Standard orientation:"]:
                    ln += 5

                    while not lines[ln].strip().startswith('----'):
                        n_ats += 1
                        ln += 1
            elif line.strip() == "SCF GIAO Magnetic shielding tensor (ppm):":
                at = 0
                while at < n_ats:
                    s_line = lines[ln].split()
                    skip = False
                    try:
                        at = int(s_line[0])
                        shifts.append(float(s_line[4]))
                        
                        element = s_line[1]
                        elements.append(element)
                    except:
                        skip = True
                        
                    if not skip:
                        try:
                            spin = self.spin_dict[element][0]
                        except:
                            spin = self._ask_spin_abundance(element)
                            self.spin_dict[element] = [spin, 1]
                        
                        spins.append(spin)
                        
                    
                    ln += 1
                    
            elif line.strip() == "Total nuclear spin-spin coupling J (Hz):":
                ln += 2
                j_table = np.zeros((n_ats, n_ats))
                i = j = at = 0
                init_j = 0
                while i < n_ats and j < n_ats:
                    at = 0
                    while at < n_ats:

                        j = init_j

                        s_line = lines[ln].split()
                        at += 1

                        try:
                            i = int(s_line[0]) - 1
                        except ValueError:
                            break

                        for j_el in s_line[1:]:
                            coupling = float(j_el.replace("D", "E"))
                            j_table[i][j] = j_table[j][i] = abs(coupling)
                            j += 1

                        if i + 1 >= n_ats:
                            ln += 1
                            break

                        ln += 1
                    ln += 1
                    try:
                        init_j = int(lines[ln].split()[0]) - 1

                    except ValueError:
                        break

            ln += 1
    
        for at in range(n_ats):
            
            pd = {}
            
            pd["Element"] = elements[at]
            pd["Shift"  ] = shifts[at]
            
            try:
                j_list = []
                for j, el in enumerate(j_table[at]):
                    j_list.append([el, spins[j], elements[j], 1])
                    
                j_list = sorted(j_list, key = lambda x: x[0])
                pd["J Coupling"] = j_list
            except NameError:
                pd["J Coupling"] = {}
                
            self.peak_dict[at] = pd
            
    def set_j_degeneracy(self):
                
        for at, pd in self.peak_dict.items():
            j_list = pd["J Coupling"]
            degeneracy_j_list = []
            
            for c, s, e, d in j_list:
                if c > self.coupling_degeneracy_threshold:
                    skip = False
                    for i, (dc, ds, de, dd) in enumerate(degeneracy_j_list):
                        if abs(c - np.average(dc)) <= self.coupling_degeneracy_threshold and e == de:
                            degeneracy_j_list[i][0].append(c)
                            degeneracy_j_list[i][3] += 1
                            skip = True
                            break
                            
                    if not skip:
                        degeneracy_j_list.append([[c], s, e, 1])
                        
                        
            degeneracy_j_list = [[np.average(dc), ds, de, dd] for dc, ds, de, dd in degeneracy_j_list]
            
            self.peak_dict[at]["J Coupling"] = degeneracy_j_list
            
    
    def _ask_spin_abundance(self, element):
        try:
            input = raw_input
        except NameError:
            pass
        
        while True:
            spin = input("Input nuclear spin for element {}:".format(element))
        
            try:
                spin = float(spin)
                if spin % 0.5 == 0 and spin >= 0:
                    break
            except:
                pass
            
            print("Spin must be positive half-integer")
        
        while True:
            abundance = input("Input abundance (0-1) for element {}:".format(element))
        
            try:
                abundance = float(abundance)
                if 0 <= abundance <= 1:
                    break
            except:
                pass
            
            print("Abundance must be between 0 and 1")
            
        return [spin, abundance]

        
class EntryPopup(object):
    def __init__(self, parent, text):
        
        self.root = tk.Toplevel(parent.root)
        self.parent = parent        
        self.value = ""
        
        self.label = _add_label(self.root, {"text": text}, gk('00'))
        self.entry = _add_entry(self.root, "", {}, gk('01'))
        self.ok_button     = _add_button(self.root, {"text": "OK"    }, gk('10'), {"<Button-1>": self._ok})
        self.cancel_button = _add_button(self.root, {"text": "Cancel"}, gk('11'), {"<Button-1>": self._cancel})
        
        self.root.protocol("WM_DELETE_WINDOW", self._cl)

    def _cl(self, *args):
        self.root.destroy()
        
    def _ok(self, *args):
        self.value = self.entry.get()
        self.root.destroy()
        
    def _cancel(self, *args):
        self.root.destroy()
        
class LoadGaussianPopup(object):
    def __init__(self, parent, parser):
        
        self.root = tk.Toplevel(parent.root)
        self.parent = parent
        self.parser = parser
        
        ell = []
        for at, pd in parser.peak_dict.items():
            element = pd['Element']
            
            if element not in ell:
                ell.append(element)
                
        
        self.element_label = _add_label(self.root, {"text": "Select Element:"}, gk('00w'))
        self.elements_om, self.element = _add_optionmenu(self.root, 'H' if 'H' in ell else ell[0], ell, {}, gk('01ew'))
        
        self.reference_label = _add_label(self.root, {"text": "Reference Shift (ppm):"}, gk('10w'))
        self.reference_entry = _add_entry(self.root, "", {}, gk('11w'))
        
        self.degeneracy_label = _add_label(self.root, {"text": "Degeneracy Threshold (Hz):"}, gk('20w'))
        self.degeneracy_entry = _add_entry(self.root, "1", {}, gk('21w'))
        
        self.decouple_label = _add_label(self.root, {"text": "Decouple Elements?"}, gk('30w'))
        self.decouple       = tk.BooleanVar(value=True)
        _add_checkbutton(self.root, True, {}, gk('31'), variable=self.decouple)
        
        self.go_button = _add_button(self.root, {"text": "Go"}, gk('40ew'), {"<Button-1>": self.go})
        
        self.root.protocol("WM_DELETE_WINDOW", self._cl)

    def _cl(self, *args):
        self.root.destroy()

    def _get_ref(self):
        
        try:
            reference = float(self.reference_entry.get())
            if reference > 0:
                return reference
        except:
            pass
        
        msgbox.showerror("Error", "Reference Shift must be a positive float")
        return None
        
    def _get_degeneracy(self):
        
        try:
            degeneracy = float(self.degeneracy_entry.get())
            if degeneracy > 0:
                return degeneracy
        except:
            pass
        
        msgbox.showerror("Error", "Degeneracy Threshold must be a positive float")
        return None
        
    def go(self, *args):
        
        reference = self._get_ref()
        if reference is None:
            return

        degeneracy_threshold = self._get_degeneracy()
        if degeneracy_threshold is None:
            return
            
        chosen_element = self.element.get()
        decouple = self.decouple.get()

        self.parser.coupling_degeneracy_threshold = degeneracy_threshold
        self.parser.set_j_degeneracy()
        
        nmr_list = []
            
        for at, pd in self.parser.peak_dict.items():
            
            j_list= pd['J Coupling']
            shift = pd['Shift']
            element = pd['Element']
            
            if element == chosen_element:
                nmr_list.append([reference - shift, [[c, s, d] for c, s, e, d in j_list if not decouple or e == element]])
                
        self.parent.gaussian_nmr_list = nmr_list
        
        self.root.destroy()

def gk(string):
    grid   = "".join([s for s in string if s.isdigit()])
    sticky = "".join([s for s in string if s in "news"])
    grid = grid.ljust(6, '0')
    r,c,rs,cs,px,py = [int(s) for s in grid]
    g = {"row": r, "column": c}
    if rs: g["rowspan"]    = rs
    if cs: g["columnspan"] = cs
    if px: g["padx"]       = px
    if py: g["pady"]       = px

    if sticky: g["sticky"]   = sticky
    
    return g
        
def _add_frame(frame_kwargs={}, grid_kwargs={}):
    """Insert a frame (box) into parent.
    With text, a labelled frame is used"""
    
    if "text" in frame_kwargs:
        frame = tk.LabelFrame(**frame_kwargs)
    else:
        frame = tk.Frame(**frame_kwargs)
        
    frame.grid(**grid_kwargs)
    return frame
    
def _add_label(frame, text_kwargs={}, grid_kwargs={}, config_kwargs={}):
    """Insert a label"""
    label = tk.Label(frame, **text_kwargs)
    label.grid(**grid_kwargs)
    label.config(**config_kwargs)
    return label
    
def _add_scale(frame, val, scale_kwargs={}, grid_kwargs={}, config_kwargs={}):
    """Insert a scrollable bar"""
    variable = tk.StringVar()
    variable.set(val)
    
    scale = tk.Scale(frame, **scale_kwargs)
    scale.set(variable.get())
    scale.grid(**grid_kwargs)
    scale.config(**config_kwargs)
    scale.grid_columnconfigure(0, weight = 1)
    return scale

def _add_button(frame, button_kwargs={}, grid_kwargs={}, bind_kwargs={}, config_kwargs={}):
    "Insert a button"""
    button = tk.Button(frame, **button_kwargs)
    button.grid(**grid_kwargs)
    for k, v in bind_kwargs.items():
        button.bind(k, v)
    button.config(bg = "blue", **config_kwargs)
    
    return button
    
def _add_entry(frame, val, entry_kwargs={}, grid_kwargs={}, config_kwargs={}, attach_func=None):
    """Add a text entry"""
    variable = tk.StringVar()
    variable.set(val)
    
    entry = tk.Entry(frame, textvariable=variable, **entry_kwargs)
    entry.bind("<FocusOut>", attach_func)
    entry.grid(**grid_kwargs)
    entry.config(**config_kwargs)
    return entry
    
def _add_optionmenu(frame, val, items, optionmenu_kwargs={}, grid_kwargs={}, config_kwargs={}):
    """Add a dropdown menu"""
    variable = tk.StringVar()
    variable.set(val)
    
    optionmenu = tk.OptionMenu(frame, variable, *items, **optionmenu_kwargs)
    optionmenu.grid(**grid_kwargs)
    optionmenu.config(**config_kwargs)
    return optionmenu, variable
    
def _add_radio(frame, val, radio_kwargs={}, grid_kwargs={}, config_kwargs={}, variable=None):
    """Add a radio button"""
    if variable is None:
        variable = tk.StringVar()
        variable.set(val)
    
    radio  = tk.Radiobutton(frame, variable=variable, **radio_kwargs)
    radio.grid(**grid_kwargs)
    radio.config(**config_kwargs)
    
def _add_checkbutton(frame, val, checkbutton_kwargs={}, grid_kwargs={}, config_kwargs={}, variable=None):
    """Add a radio button"""
    if variable is None:
        variable = tk.BooleanVar()
        variable.set(val)
    
    checkbutton  = tk.Checkbutton(frame, variable=variable, **checkbutton_kwargs)
    checkbutton.grid(**grid_kwargs)
    checkbutton.config(**config_kwargs)
    return checkbutton
    
def _add_mpl_canvas(frame, figure, grid_kwargs={}):
    
    canvas = FigureCanvasTkAgg(figure, frame)
    canvas.show()
    widget = canvas.get_tk_widget()
    widget.grid(**grid_kwargs)
    return canvas
        
if __name__ == "__main__":
    
    gui = GUI()
