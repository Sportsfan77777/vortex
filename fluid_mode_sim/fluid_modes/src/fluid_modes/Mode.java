package fluid_modes;

/**
 * creates a perturbation mode with 'm' blobs
 *
 */
public class Mode {
	
	private int mode_number;
	private Blob[] blobs;
	
	/**
	 * Perturbation Mode w/ an initial orientation
	 * @param m the mode number (must be positive integer)
	 * @param angle initial angle (in degrees) of Blob 0
	 */
	public Mode(int centerX, int centerY, int m, double angle) {
		this.mode_number = m;
		
		blobs = new Blob[m];
		double theta = (Math.PI / 180.0) * angle;
		for (int i = 0; i < m; i++) {
			double x = Math.cos(theta);
			double y = Math.sin(theta);
			blobs[i] = new Blob(centerX + x, centerY + y);
			
			theta += 2 * Math.PI / m;
		}
	}
	
	public Blob[] getBlobs() {
		return blobs;
	}
}
