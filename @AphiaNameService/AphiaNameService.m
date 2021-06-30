function obj = AphiaNameService

obj.endpoint = 'https://www.marinespecies.org/aphia.php?p=soap';
obj.wsdl = 'file:///D:/MATLAB/UVP_tools/WoRMS_Aphia.xml';

obj = class(obj,'AphiaNameService');

