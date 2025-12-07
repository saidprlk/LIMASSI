%     LIMASSI: Lightweight Image ASSImilation Framework for remite sensing 
%     to study nearshore processes

%     Copyright (C) 2025  M. Said PARLAK
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


raw_names1=dir(fullfile("xbr_images/raw_radar_04_04_2023/",'*imm*'));
sample_num=63; %sample number in each sequence

% get structure of analysis results
lms_xbr = assimilate_images_alt(raw_names1,sample_num,1);

