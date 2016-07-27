#!/usr/bin/python
import os
import sys
import math
from sets import Set
import re
import tkMessageBox
import Tkinter, Tkconstants, tkFileDialog
#import UCSD_Proteogenomics_Post_Process_Modules

class SpliceDB_GUI(Tkinter.Frame):
	
	def __init__(self, root):
		#python buildSpliceGraph.py
		#-a [rna dir]
		#-p [splice graph dir/output spl file name]
		#-f [dna dir]
		#-s [sixframe dir]
		
		self.rna_sam_or_spl_dir = '' 
		self.output_dir = 0 
		self.dna_dir = '' 
		self.sixframe_dir = '' 
		
		Tkinter.Frame.__init__(self, root)

		# options for buttons
		#button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}
		
		# define options for opening or saving a file
		self.file_opt = options = {}
		#options['defaultextension'] = '.txt'
		#options['filetypes'] = [('all files', '.*'), ('text files', '.txt')]
		#options['initialdir'] = 'C:\\'
		#options['initialfile'] = 'myfile.txt'
		#options['parent'] = root
		#options['title'] = 'Result file'
				
		# define buttons
		button_row_ind = 0
		Tkinter.Label(self, text='SAM(or .spl) file directory').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select directory', command=self.ask_dir_sam).grid(row=button_row_ind,column=1)
		
		button_row_ind += 1
		Tkinter.Label(self, text='Directory to write SpliceDB').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select directory', command=self.ask_dir_output).grid(row=button_row_ind,column=1)

		#self. = '' 		
		button_row_ind += 1
		Tkinter.Label(self, text='DNA file(.fasta) directory').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select directory', command=self.ask_dir_dna).grid(row=button_row_ind,column=1)

		#self. = '' 		
		button_row_ind += 1
		Tkinter.Label(self, text='Directory to write SixframeDB').grid(row=button_row_ind,column=0)
		Tkinter.Button(self, text='Select directory', command=self.ask_dir_sixframe).grid(row=button_row_ind,column=1)

		button_row_ind += 1
		Tkinter.Button(self, text ='Run', command = self.submit).grid(row=button_row_ind,column=0)

	def ask_dir_sam(self):
		self.rna_sam_or_spl_dir = tkFileDialog.askdirectory(title="Select directory")
		Tkinter.Label(self, text=self.rna_sam_or_spl_dir).grid(row=0,column=2)

	def ask_dir_output(self):
		self.output_dir = tkFileDialog.askdirectory(title="Select directory")
		Tkinter.Label(self, text=self.output_dir).grid(row=1,column=2)

	def ask_dir_dna(self):
		self.dna_dir = tkFileDialog.askdirectory(title="Select directory")
		Tkinter.Label(self, text=self.dna_dir).grid(row=2,column=2)

	def ask_dir_sixframe(self):
		self.sixframe_dir = tkFileDialog.askdirectory(title="Select directory")
		Tkinter.Label(self, text=self.sixframe_dir).grid(row=3,column=2)

	def submit(self):

#		print ''
#		print 'rna_sam_or_spl_dir:', self.rna_sam_or_spl_dir
#		print 'output_dir:', self.output_dir
#		print 'dna_dir:', self.dna_dir
#		print 'sixframe_dir:', self.sixframe_dir
		
		cmd_line_str = "python buildSpliceGraph " 
		cmd_line_str = cmd_line_str + " -a " + self.rna_sam_or_spl_dir
		cmd_line_str = cmd_line_str + " -p " + self.output_dir
		cmd_line_str = cmd_line_str + " -f " + self.dna_dir
		cmd_line_str = cmd_line_str + " -s " + self.sixframe_dir

		#print cmd_line_str
		#os.system(cmd_line_str)

#		SpliceDBModule = UCSD_SpliceDB_Modules.SpliceDB_Module(
#		, self.rna_sam_or_spl_dir
#		, self.output_dir
#		, self.dna_dir
#		, self.sixframe_dir)
#		SpliceDBModule.Process()
#
#		Tkinter.Label(self, text='Finished').grid(row=4,column=0)


if __name__=='__main__':

	root = Tkinter.Tk()
	root.wm_title("SpliceDB in GUI")
	SpliceDB_GUI(root).pack()
	root.mainloop()




