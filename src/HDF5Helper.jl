#=

    Module comprising of helper functions for saving time domain simulation data to HDF5 files.

=#

module HDF5Helper
using HDF5

"""
    create_group(file::HDF5.File, group_name::String)
Creates a group in an opened HDF5 file with name group_name. 
"""
function create_file_group!(file::HDF5.File, group_name::String)
    create_group(file, group_name);
end

"""
    create_dataset(file::HDF5.Group, dataset_name::String, datatype::DataType, dims::Tuple{Vararg{Int, N}}, chunk::Tuple{Vararg{Int, N}})
"""
function create_dataset!(file::HDF5.File, group_name::String, dataset_name::String, datatype::DataType, chunk_size::Int, num_time_steps::Int, N::Int)
    create_dataset(file[group_name], dataset_name, datatype, ((6,num_time_steps, N, N),(6,num_time_steps, N, N)), chunk = (6, chunk_size, N, N));
end

"""
    create_dataset!(file::HDF5.File, group_name::String, dataset_name::String, datatype::DataType, dataspace::Any, chunk::Tuple{Vararg{Int}})
Creates a dataset in an opened HDF5 file with path file[group_name/dataset_name].

Arguments:
file::HDF5.File: opened HDF5 file
group_name::String: name of the group
dataset_name::String: name of the dataset
datatype::DataType: datatype of the dataset
dataspace::Any: dataspace of the dataset. Structure is (dims0, dims_max), where dims0 is the initial dimension of the dataset and dims_max is the maximum dimension of the dataset.
chunk::Tuple{Vararg{Int}}: chunk dimensions of the data which will be saved to the dataset incrementally
"""
function create_dataset!(file::HDF5.File, group_name::String, dataset_name::String, datatype::DataType, dataspace::Any, chunk::Tuple{Vararg{Int}})
    create_dataset(file[group_name], dataset_name, datatype, dataspace, chunk = chunk);
end

function save_to_file_BDNK!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int)
    file["solution"]["du"][idx0:idx1, :] = data[1:end-1, 1, :];
    file["solution"]["deps"][idx0:idx1, :] = data[1:end-1, 2, :];
    file["solution"]["dlambda"][idx0:idx1, :] = data[1:end-1, 3, :];
    file["solution"]["du_dot"][idx0:idx1, :] = data[1:end-1, 4, :];
    file["solution"]["deps_dot"][idx0:idx1, :] = data[1:end-1, 5, :];
    file["solution"]["dlambda_dot"][idx0:idx1, :] = data[1:end-1, 6, :];
    file["solution"]["Constrained_IR1"][idx0:idx1, :] = data[1:end-1, 7, :];
    file["solution"]["Constrained_IR2"][idx0:idx1, :] = data[1:end-1, 8, :];
    file["solution"]["Wave_IR1"][idx0:idx1, :] = data[1:end-1, 9, :];
    file["solution"]["Wave_IR2"][idx0:idx1, :] = data[1:end-1, 10, :];
    file["solution"]["Wave_IR3"][idx0:idx1, :] = data[1:end-1, 11, :];
    flush(file)
end

function save_final_chunk_to_file_BDNK!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int, chunk_idx::Int, saved_time::Vector{Float64}, eq_resids::Vector{Float64})
    file["solution"]["t"] = saved_time;
    file["solution"]["eq_resids"] = eq_resids;
    file["solution"]["du"][idx0:idx1, :] = data[1:chunk_idx, 1, :];
    file["solution"]["deps"][idx0:idx1, :] = data[1:chunk_idx, 2, :];
    file["solution"]["dlambda"][idx0:idx1, :] = data[1:chunk_idx, 3, :];
    file["solution"]["du_dot"][idx0:idx1, :] = data[1:chunk_idx, 4, :];
    file["solution"]["deps_dot"][idx0:idx1, :] = data[1:chunk_idx, 5, :];
    file["solution"]["dlambda_dot"][idx0:idx1, :] = data[1:chunk_idx, 6, :];
    file["solution"]["Constrained_IR1"][idx0:idx1, :] = data[1:chunk_idx, 7, :];
    file["solution"]["Constrained_IR2"][idx0:idx1, :] = data[1:chunk_idx, 8, :];
    file["solution"]["Wave_IR1"][idx0:idx1, :] = data[1:chunk_idx, 9, :];
    file["solution"]["Wave_IR2"][idx0:idx1, :] = data[1:chunk_idx, 10, :];
    file["solution"]["Wave_IR3"][idx0:idx1, :] = data[1:chunk_idx, 11, :];
    flush(file)
end

function save_to_file_ECKART!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int)
    file["solution"]["xi"][idx0:idx1, :] = data[1:end-1, 1, :];
    file["solution"]["xi_dot"][idx0:idx1, :] = data[1:end-1, 2, :];
    file["solution"]["IR"][idx0:idx1, :] = data[1:end-1, 3, :];
    flush(file)
end

function save_final_chunk_to_file_ECKART!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int, chunk_idx::Int, saved_time::Vector{Float64})
    file["solution"]["t"] = saved_time;
    file["solution"]["xi"][idx0:idx1, :] = data[1:chunk_idx, 1, :];
    file["solution"]["xi_dot"][idx0:idx1, :] = data[1:chunk_idx, 2, :];
    file["solution"]["IR"][idx0:idx1, :] = data[1:chunk_idx, 3, :];
    flush(file)
end

module Cowling
using HDF5

"""
    create_group(file::HDF5.File, group_name::String)
Creates a group in an opened HDF5 file with name group_name. 
"""
function create_file_group!(file::HDF5.File, group_name::String)
    create_group(file, group_name);
end

"""
    create_dataset(file::HDF5.Group, dataset_name::String, datatype::DataType, dims::Tuple{Vararg{Int, N}}, chunk::Tuple{Vararg{Int, N}})
"""
function create_dataset!(file::HDF5.File, group_name::String, dataset_name::String, datatype::DataType, chunk_size::Int, num_time_steps::Int, N::Int)
    create_dataset(file[group_name], dataset_name, datatype, ((6,num_time_steps, N, N),(6,num_time_steps, N, N)), chunk = (6, chunk_size, N, N));
end

"""
    create_dataset!(file::HDF5.File, group_name::String, dataset_name::String, datatype::DataType, dataspace::Any, chunk::Tuple{Vararg{Int}})
Creates a dataset in an opened HDF5 file with path file[group_name/dataset_name].

Arguments:
file::HDF5.File: opened HDF5 file
group_name::String: name of the group
dataset_name::String: name of the dataset
datatype::DataType: datatype of the dataset
dataspace::Any: dataspace of the dataset. Structure is (dims0, dims_max), where dims0 is the initial dimension of the dataset and dims_max is the maximum dimension of the dataset.
chunk::Tuple{Vararg{Int}}: chunk dimensions of the data which will be saved to the dataset incrementally
"""
function create_dataset!(file::HDF5.File, group_name::String, dataset_name::String, datatype::DataType, dataspace::Any, chunk::Tuple{Vararg{Int}})
    create_dataset(file[group_name], dataset_name, datatype, dataspace, chunk = chunk);
end

function save_to_file_PF!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int)
    file["solution"]["xi"][idx0:idx1, :] = data[1:end-1, 1, :];
    file["solution"]["xi_dot"][idx0:idx1, :] = data[1:end-1, 2, :];
    file["solution"]["IR"][idx0:idx1, :] = data[1:end-1, 3, :];
    flush(file)
end

function save_final_chunk_to_file_PF!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int, chunk_idx::Int, saved_time::Vector{Float64})
    file["solution"]["t"] = saved_time;
    file["solution"]["xi"][idx0:idx1, :] = data[1:chunk_idx, 1, :];
    file["solution"]["xi_dot"][idx0:idx1, :] = data[1:chunk_idx, 2, :];
    file["solution"]["IR"][idx0:idx1, :] = data[1:chunk_idx, 3, :];
    flush(file)
end

function save_to_file_BDNK!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int)
    file["solution"]["du"][idx0:idx1, :] = data[1:end-1, 1, :];
    file["solution"]["deps"][idx0:idx1, :] = data[1:end-1, 2, :];
    file["solution"]["du_dr"][idx0:idx1, :] = data[1:end-1, 3, :];
    file["solution"]["deps_dr"][idx0:idx1, :] = data[1:end-1, 4, :];
    file["solution"]["du_dot"][idx0:idx1, :] = data[1:end-1, 5, :];
    file["solution"]["deps_dot"][idx0:idx1, :] = data[1:end-1, 6, :];
    file["solution"]["IR1"][idx0:idx1, :] = data[1:end-1, 7, :];
    file["solution"]["IR2"][idx0:idx1, :] = data[1:end-1, 8, :];
    flush(file)
end

function save_final_chunk_to_file_BDNK!(file::HDF5.File, data::AbstractArray, idx0::Int, idx1::Int, chunk_idx::Int, saved_time::Vector{Float64})
    file["solution"]["t"] = saved_time;
    file["solution"]["du"][idx0:idx1, :] = data[1:chunk_idx, 1, :];
    file["solution"]["deps"][idx0:idx1, :] = data[1:chunk_idx, 2, :];
    file["solution"]["du_dr"][idx0:idx1, :] = data[1:chunk_idx, 3, :];
    file["solution"]["deps_dr"][idx0:idx1, :] = data[1:chunk_idx, 4, :];
    file["solution"]["du_dot"][idx0:idx1, :] = data[1:chunk_idx, 5, :];
    file["solution"]["deps_dot"][idx0:idx1, :] = data[1:chunk_idx, 6, :];
    file["solution"]["IR1"][idx0:idx1, :] = data[1:chunk_idx, 7, :];
    file["solution"]["IR2"][idx0:idx1, :] = data[1:chunk_idx, 8, :];
    flush(file)
end

end
end