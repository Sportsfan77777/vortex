package fluid_modes;

import java.awt.Color;

/**
 * marks place in orbit where fluid element encountered perturbation
 *
 */
public class Collision extends Element {

	public Collision(double x, double y, double freq, Display d, double radius) {
		super(x, y, freq, d);
		this.radius = (int)radius;
		this.color = Color.CYAN;
	}

}
