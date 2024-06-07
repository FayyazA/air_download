from medio import read_img, save_img
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import pydicom
import scipy
import json
import math
import glob
import os
import h5py
import tqdm

from nibabel.affines import apply_affine 
import numpy.linalg as npl
import nibabel.processing
import nibabel as nib
from PIL import Image
from scipy.ndimage import zoom

#%matplotlib inline
storage_dir = "/data/psma-pet/storage/"

project_dir  = "/data/psma-pet/data/PETCT/DATASETS/DATASET_3_4_041122_165/data"
nifti_dir    = "{}/preprocessed/".format(storage_dir)
label_dir    = "/data/psma-pet/data/PETCT/DATASETS/DATASET_3_4_041122_165/data"
hdf5_dir = "{}/hdf5_4/".format(storage_dir)
fig_dir = "{}/fig/".format(storage_dir)
PT_spacing = []
CT_spacing = []
mask_spacing = []

excel_file_1 = "/data/psma-pet/data/PETCT/DATASETS/DATASET_3_4_041122_165/DATASET3_candidates__local_disease-log.xlsx"
excel_file_2 = "/data/psma-pet/data/PETCT/DATASETS/DATASET_3_4_041122_165/DATASET4_candidates__BCR_pos_metastatic__sm -log.xlsx"


excel_1 = pd.read_excel(excel_file_1,header=0)
excel_2 = pd.read_excel(excel_file_2,header=0)

useful_data_1 = excel_1.get(["Anon Acc #", "Anon Report Text"])
useful_data_2 = excel_2.get(["Anon Acc #", "Anon Report Text"])

assert os.path.isdir(project_dir)
assert os.path.isdir(label_dir)

[os.makedirs(dir_name, exist_ok=True) for dir_name in [project_dir, nifti_dir, hdf5_dir, fig_dir]]




# cohort_fname = os.path.join(project_dir, 'PSMA_PET_MSK_Pos_DATASET_1_50_studies_w_notes.csv')
# cohort_fname = os.path.join(project_dir, 'PSMA_PET_MSK_Pos_DATASET_2_50_studies_w_notes.csv')
# cohort_comp = pd.read_csv(cohort_fname, dtype={'Anon Acc #':np.int32, 'orig_acc_num': np.int32}, index_col=0)
# cohort_comp


def Write_Read_PT_CT_Nifties (PT_DICOM_dir, CT_DICOM_dir, anon_acc_num, nifti_dir):

    # PT ##############################################
    # Read DICOM array + affine
    PT_array, PT_metadata = read_img(PT_DICOM_dir, backend='itk', desired_ornt='RAI')
    global PT_spacing
    PT_spacing = PT_metadata.spacing
    
    # Save to Nifti w proper scan orientation
    PT_nii_pathname = os.path.join(nifti_dir, (anon_acc_num + '_PT.nii'))
    
    PT_nib_data = nib.Nifti1Image(PT_array, affine=PT_metadata.affine)
    nib.save(PT_nib_data, PT_nii_pathname)
    
    # CT ##############################################
    # Read DICOM array + affine
    CT_array, CT_metadata = read_img(CT_DICOM_dir, backend='itk', desired_ornt='RAI')
    global CT_spacing
    CT_spacing = CT_metadata.spacing
    
    # Save to Nifti w proper scan orientation
    CT_nii_pathname = os.path.join(nifti_dir, (anon_acc_num + '_CT.nii'))
    
    CT_nib_data = nib.Nifti1Image(CT_array, affine=CT_metadata.affine)
    nib.save(CT_nib_data, CT_nii_pathname)
    
    return (PT_array, PT_metadata.affine)

def Determine_Radius_of_Sphere (ROI):
    
    # Radius not given explicitly; determine from 2D circle
    x_coords = []
    y_coords = []
    z_coords = []

    for point in ROI['Point_mm']:
        x = float(point.split(', ')[0][1:])
        y = float(point.split(', ')[1])
        z = float(point.split(', ')[2][:-1])
    
        x_coords.append(x)
        y_coords.append(y)
        z_coords.append(z)

    radius = round((max(x_coords) - min(x_coords))/4 + \
                   (max(z_coords) - min(z_coords))/4, 2)
    
    return (radius)

def Determine_Center_of_Sphere (ROI):
    
    center_str = ROI['Center']

    x = float(center_str.split(', ')[0][1:])
    y = float(center_str.split(', ')[1])
    z = float(center_str.split(', ')[2][:-1])

    center = (x, y, z)
    
    return (center)

def Mask_Points_Inside_Spherical_ROI (PT_lesion_mask, center, radius, PT_affine, PT_array):
    
    x_min = center[0] - radius
    x_max = center[0] + radius

    y_min = center[1] - radius
    y_max = center[1] + radius

    z_min = center[2] - radius
    z_max = center[2] + radius


    for x in np.arange(x_min, x_max, 0.5):
        for y in np.arange(y_min, y_max, 0.5):
            for z in np.arange(z_min, z_max, 0.5):
            
                # Determine if point is inside or outside sphere
                inside_sphere = (x - center[0])**2 + (y - center[1])**2 + \
                                (z - center[2])**2 - radius**2 
                
                if inside_sphere <= 0:
            
                    i, j, k, _ = np.matmul(np.linalg.pinv(PT_affine), np.array([x, y, z, 1]))
            
                    x_index = int(round(i, 0))
                    y_index = int(round(j, 0))  
                    z_index = int(round(k, 0))
    
                    PT_lesion_mask[x_index, y_index, z_index] = 1
    
    return (PT_lesion_mask)

def Create_Binary_Lesion_Mask (PT_array, PT_affine, label_data):
    
    PT_lesion_mask = np.zeros(PT_array.shape)

    for image in label_data['Images']:

        image_ROIs = image['ROIs']

        if image_ROIs != []:
            print(datetime.datetime.now(), 
                  '  | ImageIndex: ', image['ImageIndex'], 
                  'Number of ROIs:', len(image['ROIs']))      

            image_index = image['ImageIndex']

            for ROI in image['ROIs']:

                radius = Determine_Radius_of_Sphere (ROI)
                center = Determine_Center_of_Sphere (ROI)

                print(datetime.datetime.now(), 
                      '  | Radius: ', radius, 
                      'Center: ', center)

                PT_lesion_mask = Mask_Points_Inside_Spherical_ROI (PT_lesion_mask, center, radius, \
                                                                   PT_affine, PT_array)
            
    return PT_lesion_mask

def Save_Label_Mask_to_Nifti (PT_lesion_mask, PT_affine, nifti_dir, anon_acc_num):

    mask_pathname = os.path.join(nifti_dir, (anon_acc_num + '_labels.nii'))

    nib_labels = nib.Nifti1Image(PT_lesion_mask, affine=PT_affine)
    nib.save(nib_labels, mask_pathname)

def preprocess_datetime (datetime):
    if '.' in datetime:
        datetime=datetime.split('.')[0]
    return (datetime)

def determine_SUV_conversion_factor(PT_metadata):
    
    weight_kg    = PT_metadata[0x0010, 0x1030].value
    height_m    = PT_metadata[0x0010, 0x1020].value
    dose         = PT_metadata[0x0054, 0x0016].value[0][0x0018, 0x1074].value
    
    inj_datetime = PT_metadata[0x0054, 0x0016].value[0][0x0018, 0x1072].value
    acq_datetime = PT_metadata[0x0008, 0x0032].value
    halflife_sec = PT_metadata[0x0054, 0x0016].value[0][0x0018, 0x1075].value
    
    inj_datetime = preprocess_datetime (inj_datetime)  
    acq_datetime = preprocess_datetime (acq_datetime)  
    
    elapsed_time = datetime.datetime.strptime(acq_datetime,'%H%M%S') - \
                    datetime.datetime.strptime(inj_datetime,'%H%M%S')
    elapsed_time_min = elapsed_time.total_seconds()/60 

    halflife_min   = halflife_sec/60
    weight_g       = weight_kg*1000 
    dose_corrected = dose*math.exp(-float(elapsed_time_min)*np.log(2)/halflife_min)

    SUV_conversion_factor = weight_g/dose_corrected
    
    return (SUV_conversion_factor, weight_kg, height_m)

def coord2indx(center, PT_affine):
    x = center[0]
    y = center[1]
    z = center[2]

    i, j, k, _ = np.matmul(np.linalg.pinv(PT_affine), np.array([x, y, z, 1]))
    x_ind = int(round(i, 0))
    y_ind = int(round(j, 0))  
    z_ind = int(round(k, 0))
    indices = (x_ind,y_ind,z_ind)
    return (indices)

def create_individual_mask (PT_array, PT_affine, label_data, PT_volume_suv, CT_volume,known_id ):
    all_lesions =[]
    all_mask1 = []
    all_mask1_info = []
    lesion_ID = 0
    for image in label_data['Images']:

            image_ROIs = image['ROIs']

            if image_ROIs != []:
                
                for ROI in image['ROIs']:
                    my_lesion = {}
                    my_lesion['ID'] = lesion_ID
                    my_lesion['radius'] = radius = Determine_Radius_of_Sphere (ROI)
                    my_lesion['center'] = center = Determine_Center_of_Sphere (ROI)
                    my_lesion['center_ijk'] = slicenum = coord2indx(center, PT_affine)
                    my_lesion['mask_spherical'] = Mask_Points_Inside_Spherical_ROI(np.zeros(PT_array.shape), center, radius, PT_affine, PT_array)
                    all_lesions.append(my_lesion)

                    mask1 = Mask_Points_Inside_Spherical_ROI(np.zeros(PT_array.shape), center, radius, PT_affine, PT_array)
                    mask1_info = [lesion_ID, *my_lesion['center_ijk']]
                    
                    all_mask1.append(mask1)
                    all_mask1_info.append(mask1_info)
                    lesion_ID += 1
                    

    # plt.subplots(nrows=len(all_lesions), ncols=4)
    # i = 1
    # for lesion in all_lesions:
    #     #CT image
    #     yindex = lesion['center_ijk'][1]
    #     CT_slice = CT_volume[:,yindex,:]
    #     plt.subplot(len(all_lesions), 4, 4*i-3)
    #     plt.imshow(CT_slice, cmap='gray',vmin=-1000, vmax=500)
    #     plt.tick_params(axis='x', labelsize=5) 
    #     plt.tick_params(axis='y', labelsize=5) 
    #     title_font = {'fontsize': 8, 'fontweight': 'bold'}
    #     plt.title('Slice (Coronal)'+' '+str(yindex),fontdict=title_font)
    #     #plt.title('Slice (Transverse)'+' '+str(266-74),fontdict=title_font)

    #     #PET

    #     PT_slice = PT_volume_suv[:,yindex,:]
    #     PT_rotated_image = np.rot90(PT_slice, k=1)
    #     plt.subplot(len(all_lesions), 4, 4*i-2)
    #     plt.imshow(PT_rotated_image, cmap='inferno',vmin=0, vmax=5)
    #     plt.tick_params(axis='x', labelsize=5) 
    #     plt.tick_params(axis='y', labelsize=5) 
    #     title_font = {'fontsize': 8, 'fontweight': 'bold'}
    #     plt.title('Slice (Coronal)'+' '+str(yindex),fontdict=title_font)
    #     #plt.title('Slice (Transverse)'+' '+str(266-74),fontdict=title_font)

    #     #ROI
    #     slice_roi = lesion['mask_spherical'][:,yindex,:]
    #     label_rotated_image = np.rot90(slice_roi, k=1)
    #     plt.subplot(len(all_lesions), 4, 4*i-1)
    #     plt.imshow(label_rotated_image, cmap='gray')
    #     plt.tick_params(axis='x', labelsize=5) 
    #     plt.tick_params(axis='y', labelsize=5) 
    #     title_font = {'fontsize': 8, 'fontweight': 'bold'}
    #     plt.title('ROI',fontdict=title_font)

    #     #Overlay
    #     plt.subplot(len(all_lesions), 4, 4*i)
    #     plt.imshow(PT_rotated_image, cmap='inferno',vmin=0, vmax=5)
    #     plt.imshow(label_rotated_image, cmap='jet',alpha = 0.5)
    #     plt.tick_params(axis='x', labelsize=5) 
    #     plt.tick_params(axis='y', labelsize=5) 
    #     title_font = {'fontsize': 8, 'fontweight': 'bold'}
    #     plt.title('Overlay',fontdict=title_font)

    #     i += 1
    # output_file = os.path.join(fig_dir,known_id) + '.jpg'
    # plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    return all_mask1, all_mask1_info

def create_mask2(mask1, mask1_info, PT_volume_suv):
    all_mask2 = []
    for lesion_num in range(mask1_info.shape[0]):
        PETimg = PT_volume_suv
        mask_bool = mask1[lesion_num].astype(bool)
        voi = PETimg[mask_bool]
        if voi.max() < 5:
            new_mask = (PETimg > voi.max()*0.42)*mask_bool
        else:
            new_mask = (PETimg > 3)*mask_bool

        all_mask2.append(new_mask)
    return all_mask2

def find_modality(subfolder):
    dicom_files = glob.glob(os.path.join(subfolder, '**/*.dcm'), recursive=True)
    file = dicom_files[0]
    ds = pydicom.filereader.dcmread(file)

    return ds.Modality



label_paths = glob.glob(label_dir + '/*/*.json', recursive=True)


for label_fname in tqdm.tqdm(label_paths):

    print("\n\nProcessing label filename {}.".format(label_fname))
    known_id = label_fname.split("/")[-2] 
        

    # For each file w labels
    f = open (label_fname, "r")
    label_data = json.loads(f.read())

    # anon_acc_num = label_fname.split('/')[-1].split('.')[0]
    anon_acc_num = label_fname.split('/')[-2]
    # acc_num      = cohort_comp['orig_acc_num'][cohort_comp['Anon Acc #'] == int(anon_acc_num)].iloc[0]

    print("Saving Nifti data for Anon Acc # {}.".format(anon_acc_num))
    try:
        report = useful_data_1.loc[useful_data_1["Anon Acc #"] == float(anon_acc_num), "Anon Report Text"].values[0]
    except:
        try:
            report = useful_data_2.loc[useful_data_2["Anon Acc #"] == float(anon_acc_num), "Anon Report Text"].values[0]
        except:
            report = "No report found"



    main_folder_dir_1 = label_fname.split(known_id)[0]
    main_folder_dir = os.path.join(main_folder_dir_1,known_id)+ '/' 
    all_items = os.listdir(main_folder_dir)   
    subfolder_paths = [os.path.join(main_folder_dir, item) for item in all_items if os.path.isdir(os.path.join(main_folder_dir, item))]

    for subfolder_path in subfolder_paths:
        modality = find_modality(subfolder_path)
        if modality == "PT":
            PT_DICOM_dir = subfolder_path  
        else:
            CT_DICOM_dir = subfolder_path

    # PT_DICOM_dir  = cohort_comp['PT_DICOM_dir'][cohort_comp['Anon Acc #'] == int(anon_acc_num)].iloc[0]
    # PT_DICOM_dir_1 = label_fname.split(known_id)[0]   
    # PT_DICOM_dir_2 = PT_DICOM_dir.split(known_id)[-1].strip('/')  
    # PT_DICOM_dir = os.path.join(PT_DICOM_dir_1,known_id, PT_DICOM_dir_2)+ '/'   

    # CT_DICOM_dir  = cohort_comp['CT_DICOM_dir'][cohort_comp['Anon Acc #'] == int(anon_acc_num)].iloc[0]
    # CT_DICOM_dir_1 = label_fname.split(known_id)[0]   
    # CT_DICOM_dir_2 = CT_DICOM_dir.split(known_id)[-1].strip('/')  
    # CT_DICOM_dir = os.path.join(CT_DICOM_dir_1,known_id, CT_DICOM_dir_2)+ '/'  

    PT_array, PT_affine = Write_Read_PT_CT_Nifties (PT_DICOM_dir, CT_DICOM_dir, anon_acc_num, nifti_dir)
    

    print("Processing labels for Anon Acc # {}.".format(anon_acc_num))

    PT_lesion_mask = Create_Binary_Lesion_Mask (PT_array, PT_affine, label_data)
    Save_Label_Mask_to_Nifti (PT_lesion_mask, PT_affine, nifti_dir, anon_acc_num)

    PT_fname = os.path.join(nifti_dir,known_id)+ '_PT.nii'   
    PT_volume = nib.load(PT_fname).get_fdata() 

    dicoms = glob.glob(PT_DICOM_dir + '/*.dcm')
    DICOM_tags = pydicom.dcmread(dicoms[0], stop_before_pixels=True)

    SUV_conversion_factor = determine_SUV_conversion_factor(DICOM_tags)[0]
    PT_volume_suv = PT_volume * SUV_conversion_factor

    CT_fname = os.path.join(nifti_dir,known_id)+ '_CT.nii'   
    CT_volume = nib.load(CT_fname).get_fdata()

    label_fname = os.path.join(nifti_dir,known_id)+ '_labels.nii'   
    label_volume = nib.load(label_fname).get_fdata()
    
    
    # plt.subplot(len(label_paths), 4, 4*i-3)
    # CT_rotated_image = np.rot90(CT_volume[:, 67, :], k=4)
    # plt.imshow(CT_rotated_image, cmap='gray',vmin=-1000, vmax=500)

    # plt.subplot(len(label_paths), 4, 4*i-2)
    # PT_rotated_image = np.rot90(PT_volume_suv[:, 67, :], k=1)
    # plt.imshow(PT_rotated_image, cmap='inferno',vmin=0, vmax=10)

    # plt.subplot(len(label_paths), 4, 4*i-1)
    # label_rotated_image = np.rot90(label_volume[:, 67, :], k=1)
    # plt.imshow(label_rotated_image, cmap='gray')

    # plt.subplot(len(label_paths), 4, 4*i)
    # plt.imshow(PT_rotated_image, cmap='inferno',vmin=0, vmax=5)
    # plt.imshow(label_rotated_image, cmap='jet',alpha = 0.5)
    # i += 1
    [mask1, mask1_info] = create_individual_mask (PT_array, PT_affine, label_data, PT_volume_suv, CT_volume,known_id)
    if len(mask1) == 0:
        mask1_array = []
        mask1_info_array = []
        mask2_array = []
    else: 
        mask1_array = np.stack(mask1, axis=0)
        mask1_info_array = np.stack(mask1_info, axis=0)
        mask2 = create_mask2(mask1, mask1_info_array, PT_volume_suv)
        mask2_array = np.stack(mask2, axis=0)
    
    weight = determine_SUV_conversion_factor(DICOM_tags)[2]
    height = determine_SUV_conversion_factor(DICOM_tags)[1]
    hdf5_fname = os.path.join(hdf5_dir,known_id)+ '.hdf5'  
    hf = h5py.File(hdf5_fname, 'w')
    hf['CT'] = CT_volume
    hf['PET'] = PT_volume_suv
    hf.attrs['new_resolution'] = [4,4,4]
    hf.attrs['new_CT_resolution'] = [4,4,4]
    hf['new_CT'] = zoom(CT_volume, (CT_spacing[0] / 4, CT_spacing[1] / 4, CT_spacing[2] / 4))
    hf['new_PET'] = zoom(PT_volume_suv, (PT_spacing[0] / 4, PT_spacing[1] / 4,PT_spacing[2] / 4))
    hf['mask0'] = label_volume
    hf['new_mask0'] = zoom(label_volume, (PT_spacing[0] / 4,PT_spacing[1] / 4, PT_spacing[2] / 4))
    #print(mask1_array.shape)
    #print(mask2_array.shape)
    hf['mask1'] = mask1_array
    #hf['new_mask1'] = zoom(mask1_array, (4 / PT_spacing[0], 4 / PT_spacing[1], 4 / PT_spacing[2]))
    hf['mask2'] = mask2_array
    #hf['new_mask2'] = zoom(mask2_array, (4 / PT_spacing[0], 4 / PT_spacing[1], 4 / PT_spacing[2]))
    hf['mask_info'] = mask1_info_array
    #hf['new_mask_info'] = zoom(mask1_info_array, (4 / PT_spacing[0], 4 / PT_spacing[1], 4 / PT_spacing[2]))
    hf.attrs['exam_id'] = known_id
    hf.attrs['suv_cf'] = SUV_conversion_factor
    hf.attrs['height_m'] = height
    hf.attrs['weight_kg'] = weight
    hf.attrs['report'] = report
    hf.attrs['resolution'] = PT_spacing
    hf.attrs['CT_spacing'] = CT_spacing
    TNM_stage = report[report.find("miTNM stage:")+13:report.find("miTNM stage:")+23]
    T_value = report.find("mi-T-stage:")
    N_value = report.find("mi-N-stage:")
    M_value = report.find("mi-M-stage:")
    T_stage = report[T_value+11:T_value+16]
    N_stage = report[N_value+11:N_value+16]
    M_stage = report[M_value+11:M_value+16]
    print(TNM_stage)
    print(T_stage)
    print(N_stage)
    print(M_stage)
    hf.attrs['TNM_stage'] = TNM_stage
    hf.attrs['T_stage'] = T_stage
    hf.attrs['N_stage'] = N_stage
    hf.attrs['M_stage'] = M_stage
    hf.close()





    


# # Document the location of the Nifties

# PT_NIfTI_fnames=[]
# CT_NIfTI_fnames=[]

# for index, row in cohort_comp.iterrows():
    
#     print("Processing Anon Acc # ", row['Anon Acc #'])
    
#     if row['JSON'] == 1:
#         PT_NIfTI_fname = os.path.join(nifti_dir, (str(row['Anon Acc #']) + "_PT.nii"))
#         CT_NIfTI_fname = os.path.join(nifti_dir, (str(row['Anon Acc #']) + "_CT.nii"))
#     else:
#         PT_NIfTI_fname = None
#         CT_NIfTI_fname = None
    
#     PT_NIfTI_fnames.append(PT_NIfTI_fname)
#     CT_NIfTI_fnames.append(CT_NIfTI_fname)
        
        
# cohort_comp["PT_NIfTI_fname"] = PT_NIfTI_fnames
# cohort_comp["CT_NIfTI_fname"] = CT_NIfTI_fnames

