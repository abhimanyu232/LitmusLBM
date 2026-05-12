import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
import imageio.v3 as imageio

def read_data_file(filename):
		"""Read data from a text file."""
		return np.loadtxt(filename)

def plot_velocity_field(velocity_data, timestep, output_dir='plots'):
		"""Create a contour plot of the velocity magnitude field."""
		fig = plt.figure(figsize=(10, 8), dpi=96)
		# Create a contour plot
		# levels = np.linspace(velocity_data.min(), velocity_data.max(), 50)
		levels = np.linspace(0.0, 0.12, 200)

		contour = plt.contourf(velocity_data, levels=levels, cmap='viridis')

		# Add a colorbar
		plt.colorbar(contour, label='Velocity Magnitude')

		# Add title and labels
		plt.title(f'Velocity Magnitude Field (t = {timestep})')
		plt.xlabel('x')
		plt.ylabel('y')


		# Create output directory if it doesn't exist
		os.makedirs(output_dir, exist_ok=True)

		# fig.set_size_inches(10, 8)
		# plt.gcf().set_size_inches(10, 10)

		# Save the plot
		fig.tight_layout()
		plt.savefig(os.path.join(output_dir, f'velocity_field_{timestep}.png'))
		plt.close()

def plot_vorticity_field(vorticity_data, timestep, output_dir='plots'):
		"""Create a contour plot of the vorticity field."""
		fig = plt.figure(figsize=(10, 8), dpi=96)

		# Create a contour plot with symmetric colormap for vorticity
		vmax = max(abs(vorticity_data.min()), abs(vorticity_data.max()))
		# levels = np.linspace(-vmax, vmax, 50)
		levels = np.linspace(-0.004, 0.004, 200)

		contour = plt.contourf(vorticity_data, levels=levels, cmap='RdBu_r')

		# Add a colorbar
		plt.colorbar(contour, label='Vorticity')

		# Add title and labels
		plt.title(f'Vorticity Field (t = {timestep})')
		plt.xlabel('x')
		plt.ylabel('y')

		# Create output directory if it doesn't exist
		os.makedirs(output_dir, exist_ok=True)

		# Save the plot
		fig.tight_layout()
		plt.savefig(os.path.join(output_dir, f'vorticity_field_{timestep}.png'))
		plt.close()

def plot_combined_fields(velocity_data, vorticity_data, timestep, output_dir='plots'):
		"""Create a combined plot showing both velocity and vorticity."""
		fig, (ax1, ax2) = plt.subplots(1, 2)
		fig.set_size_inches(20, 8)

		# Velocity plot
		# vel_levels = np.linspace(velocity_data.min(), velocity_data.max(), 50)
		vel_levels = np.linspace(0.0, 0.12, 200)

		vel_contour = ax1.contourf(velocity_data, levels=vel_levels, cmap='viridis')
		ax1.set_title(f'Velocity Magnitude (t = {timestep})')
		ax1.set_xlabel('x')
		ax1.set_ylabel('y')
		plt.colorbar(vel_contour, ax=ax1, label='Velocity Magnitude')

		# Vorticity plot
		vmax = max(abs(vorticity_data.min()), abs(vorticity_data.max()))
		# vort_levels = np.linspace(-vmax, vmax, 50)
		vort_levels = np.linspace(-0.004, 0.004, 200)

		vort_contour = ax2.contourf(vorticity_data, levels=vort_levels, cmap='RdBu_r')
		ax2.set_title(f'Vorticity (t = {timestep})')
		ax2.set_xlabel('x')
		ax2.set_ylabel('y')
		plt.colorbar(vort_contour, ax=ax2, label='Vorticity')

		# Create output directory if it doesn't exist
		os.makedirs(output_dir, exist_ok=True)

		# Save the combined plot
		fig.tight_layout()
		plt.savefig(os.path.join(output_dir, f'combined_fields_{timestep}.png'), dpi=96)
		plt.close()

def _numeric_key_from_filename(fname):
		base = os.path.basename(fname)
		m = re.search(r'_(\d+)\.png$', base)
		return int(m.group(1)) if m else 0

def create_mp4_from_images(image_glob, output_path, fps=6):
		files = sorted(glob.glob(image_glob), key=_numeric_key_from_filename)
		if not files:
				print(f"No images found for pattern: {image_glob}. Skipping {output_path}.")
				return
		print(f"Creating animation: {output_path} from {len(files)} frames")
		os.makedirs(os.path.dirname(output_path), exist_ok=True)
		frames = [imageio.imread(f) for f in files]
		try:
				imageio.imwrite(output_path, frames, plugin='FFMPEG', fps=fps)
		except Exception as e:
				print(f"Failed to write {output_path} with FFMPEG plugin: {e}.\n"
							f"Hint: install imageio-ffmpeg (pip install imageio-ffmpeg).")

def create_all_animations(output_dir='plots', fps=6):
		create_mp4_from_images(os.path.join(output_dir, 'velocity_field_*.png'), os.path.join(output_dir, 'velocity_animation.mp4'), fps=fps)
		create_mp4_from_images(os.path.join(output_dir, 'vorticity_field_*.png'), os.path.join(output_dir, 'vorticity_animation.mp4'), fps=fps)
		create_mp4_from_images(os.path.join(output_dir, 'combined_fields_*.png'), os.path.join(output_dir, 'combined_animation.mp4'), fps=fps)

def main():
		# Get all velocity and vorticity data files from result_fields directory
		velocity_files = sorted(glob.glob('result_fields/velocity_*.txt'))
		vorticity_files = sorted(glob.glob('result_fields/vorticity_*.txt'))

		if not velocity_files:
				print("No velocity data files found in result_fields directory!")
				return

		if not vorticity_files:
				print("No vorticity data files found in result_fields directory!")
				return

		print(f"Found {len(velocity_files)} velocity field files")
		print(f"Found {len(vorticity_files)} vorticity field files")

		# Process each timestep
		for vel_file, vort_file in zip(velocity_files, vorticity_files):
				# Extract timestep from filename
				vel_filename = os.path.basename(vel_file)
				timestep = int(vel_filename.split('_')[1].split('.')[0])

				# Read the data
				print(f"Processing timestep {timestep}...")
				velocity_data = read_data_file(vel_file)
				vorticity_data = read_data_file(vort_file)

				# Create plots
				plot_velocity_field(velocity_data, timestep)
				plot_vorticity_field(vorticity_data, timestep)	
				plot_combined_fields(velocity_data, vorticity_data, timestep)

		# for vort_file in vorticity_files:
		# 	vort_filename = os.path.basename(vort_file)
		# 	timestep = int(vort_filename.split('_')[1].split('.')[0])
		# 	print(f"Processing timestep {timestep}...")
   
		# 	vorticity_data = read_data_file(vort_file)
		# 	plot_vorticity_field(vorticity_data, timestep)

   

		# print("\nVisualization complete! Plots have been saved in the 'plots' directory.")
		print("Generated plots:")
		print("  - velocity_field_*.png: Velocity magnitude contours")
		print("  - vorticity_field_*.png: Vorticity contours")
		print("  - combined_fields_*.png: Side-by-side comparison")

		# Create MP4 animations from generated PNGs
		create_all_animations(output_dir='plots', fps=12)
		print("\nAnimations saved in 'plots':")
		print("  - velocity_animation.mp4")
		print("  - vorticity_animation.mp4")
		print("  - combined_animation.mp4")


if __name__ == "__main__":
		main()