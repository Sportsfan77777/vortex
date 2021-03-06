package fluid_modes;

import java.awt.Color;

/**
 * Blob is an Element that is a part of a Mode
 * For Mode m, there are m Blobs
 * The default appearance of a Blob is different than that of a generic Element
 */
public class Blob extends Element {
	
	public Blob(double x, double y, double freq, Display d) {
		super(x, y, freq, d);
		this.radius = 15;
		this.color = Color.BLUE;
	}
}
