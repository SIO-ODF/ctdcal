import sys
import os
import shutil
import csv
import json
import random
import pandas as pd


bower_components = os.path.join(os.path.dirname(__file__),'webviewer_Files/bower_components')
custom_js = os.path.join(os.path.dirname(__file__),'webviewer_Files/js')
index_file = os.path.join(os.path.dirname(__file__),'webviewer_Files/index.html')

dataDir = 'data'

subSampleRate = 5  #1/subSampleRate reduction in filesize

DEBUG = False

def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class WebViewer():
	"""docstring for WebView"""
	def __init__(self, parentDir='./'):
		self.parentDir = parentDir

	def _setDebug(self):
		global DEBUG
		DEBUG = True

	def _setWebviewerFolder(self, dirName):
		self.parentDir = dirName

	def _buildScaffolding(self):
		if os.path.isdir(self.parentDir):
			try:
				#debugPrint('Copy', index_file, 'to', os.path.join(self.parentDir'))
				shutil.copy(index_file,self.parentDir)
			except FileExistsError as exception:
				pass
			except PermissionError as exception:
				debugPrint('ERROR: do not have write permission at destination directory')
				return False

			try:
				#debugPrint('Copy', bower_components, 'to', os.path.join(self.parentDir, 'bower_components'))
				shutil.copytree(bower_components,os.path.join(self.parentDir, 'bower_components'))
			except FileExistsError as exception:
				pass
			except PermissionError as exception:
				debugPrint('ERROR: do not have write permission at destination directory')
				return False

			try:
				#debugPrint('Copy', bower_components, 'to', os.path.join(self.parentDir, 'bower_components'))
				shutil.copytree(custom_js,os.path.join(self.parentDir, 'js'))
			except FileExistsError as exception:
				pass
			except PermissionError as exception:
				debugPrint('ERROR: do not have write permission at destination directory')
				return False

			try:
				#debugPrint('Adding', os.path.join(self.parentDir,'data'), 'directory')
				os.mkdir(os.path.join(self.parentDir,dataDir), 0o755)
			except FileExistsError as exception:
				pass
			except PermissionError as exception:
				debugPrint('ERROR: do not have write permission at destination directory')
				return False

		else:
			errPrint('Parent Directory:', self.parentDir, 'does not exists')
			return False

		return True

	def _buildData(self, datafile):
		with open(datafile, newline='') as csvfile:
			dataReader = csv.reader(csvfile)
			raw_short_name = next(dataReader)
			raw_units = next(dataReader)
			raw_dataType = next(dataReader)

		debugPrint('raw_short_name:', raw_short_name)
		debugPrint('raw_units:', raw_units)
		debugPrint('raw_dataType:', raw_dataType)

		proc_idx = []
		proc_short_name = []
		proc_units = []
		proc_dataType = []

		for idx, val in enumerate(raw_dataType):
			if val == 'float64':
				proc_idx.append(idx)
				proc_short_name.append(raw_short_name[idx])
				proc_units.append(raw_units[idx])
				proc_dataType.append(raw_dataType[idx])

		#debugPrint('proc_idx:',proc_idx)
		#debugPrint('proc_short_name:',proc_short_name)
		#debugPrint('proc_units:',proc_units)
		#debugPrint('proc_dataType:',proc_dataType)

		df = pd.read_csv(datafile, skiprows=[1,2], usecols=proc_idx)
		#debugPrint(df.head())

		#debugPrint(df.iloc[:,0].head())

		output = {
			'castName':datafile,
			'visualizerData':[],
			'stats':[]
		}

		total_rows = len(df.index)
		#debugPrint('total_rows:', total_rows)
		
		# The row indices to skip - make sure 0 is not included to keep the header!
		skip_idx = [x for x in range(3, total_rows) if x % subSampleRate == 0]
		#debugPrint(skip_idx)

		for idx, val in enumerate(proc_short_name):
			stat = {'statName': proc_short_name[idx] + ' Bounds','statUnit': proc_units[idx], 'statType':'bounds', 'statData':[round(df.iloc[:,idx].min(),3), round(df.iloc[:,idx].max(),3)]}
			output['stats'].append(stat)

			data = {'data': df.iloc[skip_idx,idx].tolist(), 'unit':proc_units[idx], 'label':proc_short_name[idx]}
			output['visualizerData'].append(data)

		return json.dumps(output)

	def _buildDataFromDF(self, dataframe, castName):

		output = {
			'castName':castName,
			'visualizerData':[],
			'stats':[]
		}

		total_rows = len(dataframe.index)
		#debugPrint('total_rows:', total_rows)

		#get name/units
		headers = dataframe.columns.values.tolist()

		proc_units = []
		proc_short_name = []
		for header in headers:
			header_array = header.split('_')
			proc_units.append(header_array.pop())
			proc_short_name.append('_'.join(header_array))

		
		# The row indices to skip - make sure 0 is not included to keep the header!
		skip_idx = [x for x in range(3, total_rows) if x % subSampleRate == 0]
		#debugPrint(skip_idx)

		for idx, val in enumerate(proc_short_name):
			if proc_short_name != 'index':
				stat = {'statName': proc_short_name[idx] + ' Bounds','statUnit': proc_units[idx], 'statType':'bounds', 'statData':[round(dataframe.iloc[:,idx].min(),3), round(dataframe.iloc[:,idx].max(),3)]}
				output['stats'].append(stat)

				data = {'data': dataframe.iloc[skip_idx,idx].tolist(), 'unit':proc_units[idx], 'label':proc_short_name[idx]}
				output['visualizerData'].append(data)

		return json.dumps(output)

	def _saveData(self, output):
		dataFilename = 'data.json'
		dataFilePath = os.path.join(self.parentDir, dataDir, dataFilename)
		with open(dataFilePath, 'w') as f:
			f.write(output)